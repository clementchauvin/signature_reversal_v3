import os
import duckdb
import pandas as pd
from classe_phenotype import PhenotypeDefinition


class CohortExtractor:
    """Classe responsable de l'extraction des cas et des témoins éligibles

    à partir des fichiers OMOP CDM via DuckDB, incluant le comptage
    de l'épaisseur du dossier (nombre d'événements).
    """

    def __init__(self, parquet_path: str, con=None):
        """Initialise l'extracteur et connecte DuckDB au dossier contenant les

        tables OMOP.

        :param parquet_path: Dossier contenant les fichiers .parquet ou .csv.zst.
        :param con: Connexion DuckDB existante (optionnelle).
        """
        self.parquet_path = parquet_path
        # Utilisation d'une base DuckDB en mémoire pour la vitesse de calcul
        self.con = con or duckdb.connect(":memory:")
        # Création des vues virtuelles sur les tables OMOP
        self._init_omop_views()

    def _init_omop_views(self):
        """Crée des vues virtuelles DuckDB pointant vers les fichiers OMOP.

        Cela évite de charger les fichiers complets en mémoire RAM Python.
        """
        tables = [
            "person",
            "condition_occurrence",
            "drug_exposure",
            "measurement",
            "procedure_occurrence",
            "observation_period",
            "concept",
        ]
        for table in tables:
            p_parquet = os.path.join(self.parquet_path, f"{table}.parquet")
            p_csv = os.path.join(self.parquet_path, f"{table}.csv.zst")

            # Détection automatique du format disponible (.parquet ou .csv.zst)
            if os.path.exists(p_parquet):
                self.con.execute(
                    f"CREATE VIEW IF NOT EXISTS {table} AS SELECT * FROM read_parquet('{p_parquet}')"
                )
            elif os.path.exists(p_csv):
                self.con.execute(
                    f"CREATE VIEW IF NOT EXISTS {table} AS SELECT * FROM read_csv_auto('{p_csv}')"
                )

    def _build_signal_subqueries(
        self, phenotype: PhenotypeDefinition
    ) -> list[str]:
        """Méthode privée qui génère dynamiquement les sous-requêtes SQL UNION

        selon les critères définis dans le phénotype clinique.
        """
        subqueries = []

        # 1. Diagnostics (condition_occurrence)
        if phenotype.disease_concepts:
            concepts_str = ",".join(map(str, phenotype.disease_concepts))
            subqueries.append(f"""
                SELECT person_id, condition_start_date AS event_date 
                FROM condition_occurrence 
                WHERE condition_concept_id IN ({concepts_str})
            """)

        # 2. Prescriptions / Médicaments (drug_exposure)
        if phenotype.drug_concepts:
            concepts_str = ",".join(map(str, phenotype.drug_concepts))
            subqueries.append(f"""
                SELECT person_id, drug_exposure_start_date AS event_date 
                FROM drug_exposure 
                WHERE drug_concept_id IN ({concepts_str})
            """)

        # 3. Biologie / Labos (measurement) avec seuils numériques
        if phenotype.lab_conditions:
            for lab in phenotype.lab_conditions:
                subqueries.append(f"""
                    SELECT person_id, measurement_date AS event_date 
                    FROM measurement 
                    WHERE measurement_concept_id = {lab['concept_id']} 
                      AND value_as_number >= {lab['min_value']}
                """)

        # 4. Actes & Chirurgies (procedure_occurrence)
        if phenotype.procedure_concepts:
            concepts_str = ",".join(map(str, phenotype.procedure_concepts))
            subqueries.append(f"""
                SELECT person_id, procedure_date AS event_date 
                FROM procedure_occurrence 
                WHERE procedure_concept_id IN ({concepts_str})
            """)

        return subqueries

    def _build_all_events_view(self):
        """Crée une vue unifiée de TOUS les événements médicaux (toutes tables

        confondues) pour permettre le comptage rapide du nombre de lignes du
        dossier.
        """
        self.con.execute("""
            CREATE OR REPLACE VIEW all_medical_events AS
            SELECT person_id, condition_start_date AS event_date FROM condition_occurrence
            UNION ALL
            SELECT person_id, drug_exposure_start_date AS event_date FROM drug_exposure
            UNION ALL
            SELECT person_id, measurement_date AS event_date FROM measurement
            UNION ALL
            SELECT person_id, procedure_date AS event_date FROM procedure_occurrence
        """)

    def extract_cases(
        self, phenotype: PhenotypeDefinition, lookback_days: int = 365
    ) -> pd.DataFrame:
        """Extrait les cas atteints de la maladie, calcule leur t0 précoce,

        vérifie la fenêtre de suivi minimale, et compte le nombre de
        lignes/événements dans leur dossier avant t0.
        """
        # Générer les sous-requêtes spécifiques au phénotype
        subqueries = self._build_signal_subqueries(phenotype)
        if not subqueries:
            raise ValueError(
                f"Aucun critère défini dans le phénotype {phenotype.name}"
            )

        union_sql = "\nUNION ALL\n".join(subqueries)

        # Créer la vue unifiée pour le comptage du volume de données
        self._build_all_events_view()

        # Requête SQL DuckDB principale pour les CAS
        query = f"""
        WITH disease_signals AS (
            {union_sql}
        ),
        -- 1. Détermination de la date t0 la plus précoce (MIN)
        first_events AS (
            SELECT 
                person_id AS patient_id,
                MIN(event_date) AS prediction_time
            FROM disease_signals
            GROUP BY person_id
        ),
        -- 2. Filtrage des cas selon la présence d'un historique minimal (lookback_days)
        valid_cases AS (
            SELECT 
                fe.patient_id,
                fe.prediction_time,
                p.gender_concept_id AS gender,
                p.birth_datetime,
                op.observation_period_start_date,
                op.observation_period_end_date
            FROM first_events fe
            JOIN person p ON fe.patient_id = p.person_id
            JOIN observation_period op ON fe.patient_id = op.person_id
            WHERE op.observation_period_start_date <= (fe.prediction_time - INTERVAL '{lookback_days}' DAY)
              AND fe.prediction_time <= op.observation_period_end_date
        )
        -- 3. Comptage du nombre total de lignes/événements avant t0 pour chaque cas (Épaisseur du dossier)
        SELECT 
            vc.patient_id,
            vc.prediction_time,
            vc.gender,
            vc.birth_datetime,
            vc.observation_period_start_date,
            vc.observation_period_end_date,
            COUNT(e.event_date) AS event_count -- ÉPAISSEUR DU DOSSIER EN NOMBRE DE LIGNES
        FROM valid_cases vc
        LEFT JOIN all_medical_events e 
               ON vc.patient_id = e.person_id 
              AND e.event_date <= vc.prediction_time
        GROUP BY 
            vc.patient_id, 
            vc.prediction_time, 
            vc.gender, 
            vc.birth_datetime, 
            vc.observation_period_start_date, 
            vc.observation_period_end_date
        """

        return self.con.execute(query).df()

    def extract_eligible_controls(
        self, cases_df: pd.DataFrame, phenotype: PhenotypeDefinition
    ) -> pd.DataFrame:
        """Extrait l'ensemble des témoins potentiels en excluant tous les cas

        et tous les patients présentant un signal clinique de la maladie.
        """
        # Enregistrement du DataFrame des cas dans DuckDB pour l'exclusion
        self.con.register(
            "cases_to_exclude",
            pd.DataFrame({"person_id": cases_df["patient_id"]}),
        )

        subqueries = self._build_signal_subqueries(phenotype)
        exclusion_union = (
            "\nUNION DISTINCT\n".join(subqueries)
            if subqueries
            else "SELECT -1 AS person_id, '1900-01-01'::DATE AS event_date"
        )

        # Extraction des témoins sains avec leurs métadonnées
        query = f"""
        WITH all_excluded_patients AS (
            SELECT person_id FROM cases_to_exclude
            UNION DISTINCT
            SELECT person_id FROM ({exclusion_union})
        )
        SELECT 
            p.person_id AS patient_id,
            p.gender_concept_id AS gender,
            p.birth_datetime,
            op.observation_period_start_date,
            op.observation_period_end_date
        FROM person p
        JOIN observation_period op ON p.person_id = op.person_id
        WHERE p.person_id NOT IN (SELECT person_id FROM all_excluded_patients)
        """

        controls_df = self.con.execute(query).df()

        # Enregistrement temporaire des témoins dans DuckDB pour pré-calculer
        # le total d'événements global (qui servira lors de l'alignement sur le t0 du cas)
        self.con.register("temp_eligible_controls", controls_df)

        # On rajoute la vue globale si pas encore faite
        self._build_all_events_view()

        # Comptage global des lignes des témoins
        controls_with_counts_query = """
        SELECT 
            c.*,
            COUNT(e.event_date) AS total_lifetime_events -- Nombre total de lignes dans la vie du témoin
        FROM temp_eligible_controls c
        LEFT JOIN all_medical_events e ON c.patient_id = e.person_id
        GROUP BY 
            c.patient_id, c.gender, c.birth_datetime, 
            c.observation_period_start_date, c.observation_period_end_date
        """

        return self.con.execute(controls_with_counts_query).df()