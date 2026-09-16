import os
import sys
from datetime import datetime

# Import des classes du pipeline
from classe_phenotype import PhenotypeDefinition
from classe_extractor import CohortExtractor
from classe_matcher import CohortMatcher
from classe_pipeline import RABITPipeline


def main():
    # =========================================================================
    # 0. CONFIGURATION DES CHEMINS ET PARAMS
    # =========================================================================
    DATA_DIR = "./data/omop_parquet"      # Dossier contenant les fichiers .parquet OMOP CDM
    WORK_DIR = "./outputs/run_t2d"        # Dossier de travail pour les sorties intermédiaires
    FINAL_CSV = os.path.join(WORK_DIR, "delta_p_statistics.csv")
    
    RABIT_SCRIPT = "./rabit/run_rabit.py" # Chemin vers le script d'inférence RABIT
    DATA_SOURCE = "shc"                   # Source de données configurée pour RABIT
    GPU_ID = 0                            # Identifiant du GPU à utiliser
    
    K_NEIGHBORS = 5                       # Ratio 1:k (5 témoins par cas)
    LOOKBACK_DAYS = 365                   # Fenêtre de suivi minimale avant t0 (1 an)
    MAX_AGE_DIFF = 2.0                    # Écart d'âge maximal toléré (2 ans)
    MAX_EVENT_RATIO = 2.0                 # Tolérance sur le volume de lignes (entre 50% et 200%)

    os.makedirs(WORK_DIR, exist_ok=True)
    print("=" * 80)
    print(f"[{datetime.now().strftime('%H:%M:%S')}] DÉMARRAGE DU PIPELINE DE RECHERCHE PROTEOMIQUE")
    print("=" * 80)

    # =========================================================================
    # 1. DÉFINITION DU PHÉNOTYPE CLINIQUE (Exemple : Diabète de Type 2 - T2D)
    # =========================================================================
    print("\n[Étape 1/4] Configuration du phénotype clinique...")
    
    t2d_phenotype = PhenotypeDefinition(
        name="Type_2_Diabetes",
        # Concepts OMOP condition (ex: Diabète de Type 2)
        disease_concepts=[201826, 443767, 4193704],
        # Concepts OMOP médicaments (ex: Metformine, Insuline, Sulfamides)
        drug_concepts=[1503297, 1502809, 1502855],
        # Conditions biologie (ex: HbA1c >= 6.5% / Concept 3004410)
        lab_conditions=[
            {"concept_id": 3004410, "min_value": 6.5}
        ],
        # Prestations / Actes (optionnel)
        procedure_concepts=[]
    )
    print(f"-> Phénotype configuré : {t2d_phenotype.name}")

    # =========================================================================
    # 2. EXTRACTION DE LA COHORTE VIA DUCKDB (CohortExtractor)
    # =========================================================================
    print("\n[Étape 2/4] Extraction des cas et témoins éligibles via DuckDB...")
    
    extractor = CohortExtractor(parquet_path=DATA_DIR)
    
    # A. Extraction des cas avec t0 précoce et comptage de l'épaisseur du dossier (lignes)
    cases_df = extractor.extract_cases(
        phenotype=t2d_phenotype, 
        lookback_days=LOOKBACK_DAYS
    )
    print(f"-> Cas extraits : {len(cases_df)} patients répondant aux critères.")

    if cases_df.empty:
        print("[ERREUR] Aucun cas extrait. Arrêt du pipeline.")
        sys.exit(1)

    # B. Extraction de la réserve de témoins sains éligibles
    controls_df = extractor.extract_eligible_controls(
        cases_df=cases_df, 
        phenotype=t2d_phenotype
    )
    print(f"-> Témoins éligibles identifiés : {len(controls_df)} patients sains.")

    # =========================================================================
    # 3. APPARIEMENT 1:K STRUCTURÉ (CohortMatcher)
    # =========================================================================
    print("\n[Étape 3/4] Exécution de l'appariement 1:k sans remplacement...")
    
    matcher = CohortMatcher(
        k_neighbors=K_NEIGHBORS,
        max_age_diff_years=MAX_AGE_DIFF,
        max_event_count_ratio=MAX_EVENT_RATIO
    )
    
    matched_cohort_df = matcher.match(cases_df=cases_df, controls_df=controls_df)
    
    n_matched_cases = matched_cohort_df[matched_cohort_df["is_case"] == 1]["patient_id"].nunique()
    n_matched_controls = matched_cohort_df[matched_cohort_df["is_case"] == 0]["patient_id"].nunique()
    
    print(f"-> Appariement réussi !")
    print(f"   - Cas retenus : {n_matched_cases}")
    print(f"   - Témoins retenus : {n_matched_controls}")
    print(f"   - Nombre total de lignes d'inférence : {len(matched_cohort_df)}")

    # =========================================================================
    # 4. INFÉRENCE GPU RABIT & CALCUL DU DELTA P (RABITPipeline)
    # =========================================================================
    print("\n[Étape 4/4] Inférence GPU RABIT et calcul des statistiques protéomiques...")
    
    pipeline = RABITPipeline(
        rabit_script_path=RABIT_SCRIPT,
        data_source=DATA_SOURCE,
        gpu_id=GPU_ID,
        batch_size=256
    )
    
    # Exécution complète : Export CSV -> Inférence GPU -> Calcul Delta P -> Sauvegarde CSV
    stats_df = pipeline.run_full_pipeline(
        cohort_df=matched_cohort_df,
        work_dir=WORK_DIR,
        output_stats_csv=FINAL_CSV
    )

    # =========================================================================
    # RESUME DU RESULTAT
    # =========================================================================
    print("\n" + "=" * 80)
    print(f"[{datetime.now().strftime('%H:%M:%S')}] PIPELINE EXÉCUTÉ AVEC SUCCÈS")
    print("=" * 80)
    print(f"Statistiques générées pour {len(stats_df)} protéines.")
    print("Aperçu des 5 premières protéines :")
    print(stats_df.head())
    print(f"\nLivrable final disponible sous : {FINAL_CSV}")


if __name__ == "__main__":
    main()