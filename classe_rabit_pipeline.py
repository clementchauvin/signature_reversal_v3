import os
import subprocess
import numpy as np
import pandas as pd


class RABITPipeline:
    """Classe responsable de l'orchestration de l'inférence RABIT sur GPU

    et du calcul vectorisé des statistiques protéomiques (Delta P).
    """

    def __init__(
        self,
        rabit_script_path: str,
        data_source: str = "shc",
        gpu_id: int = 0,
        batch_size: int = 256,
    ):
        """Initialise la configuration d'exécution de RABIT.

        :param rabit_script_path: Chemin vers le script d'inférence de RABIT
            (ex: 'run_rabit.py').
        :param data_source: Source de données configurée dans RABIT (ex: 'shc').
        :param gpu_id: Identifiant du GPU à utiliser (ex: 0).
        :param batch_size: Taille de batch pour l'inférence GPU.
        """
        self.rabit_script_path = rabit_script_path
        self.data_source = data_source
        self.gpu_id = gpu_id
        self.batch_size = batch_size

    def export_prediction_times(
        self, cohort_df: pd.DataFrame, output_csv_path: str
    ):
        """Exporte le DataFrame d'appariement au format CSV attendu par RABIT

        (colonnes: patient_id, prediction_time).
        """
        # RABIT a besoin uniquement de l'ID patient et de l'index temporel t0
        rabit_input = cohort_df[["patient_id", "prediction_time"]].copy()

        # Formatage strict de la date au format YYYY-MM-DD
        rabit_input["prediction_time"] = pd.to_datetime(
            rabit_input["prediction_time"]
        ).dt.strftime("%Y-%m-%d")

        # Sauvegarde du CSV sans index
        rabit_input.to_csv(output_csv_path, index=False)
        print(
            f"[RABITPipeline] Fichier d'entrée exporté : {output_csv_path} ({len(rabit_input)} lignes)"
        )

    def run_inference(
        self, pred_times_csv: str, output_parquet_dir: str
    ) -> str:
        """Lance l'inférence RABIT sur GPU via un sous-processus système.

        :param pred_times_csv: Chemin vers le fichier CSV généré par
            export_prediction_times.
        :param output_parquet_dir: Dossier où RABIT doit enregistrer les
            protéomes synthétiques.
        :return: Chemin du fichier Parquet de sortie généré par RABIT.
        """
        os.makedirs(output_parquet_dir, exist_ok=True)

        # Construction de la commande shell pour exécuter RABIT sur le GPU cible
        cmd = [
            "python",
            self.rabit_script_path,
            "--pred_times",
            pred_times_csv,
            "--data_source",
            self.data_source,
            "--out_dir",
            output_parquet_dir,
            "--batch_size",
            str(self.batch_size),
        ]

        # Définition de la variable d'environnement pour sélectionner le GPU
        env = os.environ.copy()
        env["CUDA_VISIBLE_DEVICES"] = str(self.gpu_id)

        print(
            f"[RABITPipeline] Lancement de RABIT sur GPU:{self.gpu_id}..."
        )

        # Exécution de la commande et capture des erreurs éventuelles
        result = subprocess.run(
            cmd, env=env, capture_output=True, text=True, check=True
        )

        print("[RABITPipeline] Inférence terminée avec succès.")

        # Chemin du fichier Parquet généré par RABIT
        output_parquet_file = os.path.join(
            output_parquet_dir, "predictions.parquet"
        )
        return output_parquet_file

    def compute_protein_statistics(
        self, synthetic_proteome_parquet: str, cohort_df: pd.DataFrame
    ) -> pd.DataFrame:
        """Calcule de façon vectorisée les statistiques de différence

        protéique (Delta P) entre les Cas et les Témoins appariés.

        :param synthetic_proteome_parquet: Chemin du fichier Parquet de sortie
            de RABIT.
        :param cohort_df: DataFrame d'appariement contenant match_group_id,
            patient_id, is_case.
        :return: DataFrame Pandas indexé par les 2 923 protéines (colonnes:
            mean, var, std, sem).
        """
        print(
            "[RABITPipeline] Chargement des prédictions et calcul du Delta P vectorisé..."
        )

        # 1. Lecture du fichier Parquet des protéomes synthétiques
        predictions_df = pd.read_parquet(synthetic_proteome_parquet)

        # 2. Fusion avec la cohorte pour récupérer match_group_id et is_case
        df = pd.merge(
            predictions_df,
            cohort_df[["patient_id", "match_group_id", "is_case"]],
            on="patient_id",
            how="inner",
        )

        # Identification automatique de toutes les colonnes de protéines (ex: P01023, P02751...)
        metadata_cols = {
            "patient_id",
            "prediction_time",
            "match_group_id",
            "is_case",
        }
        protein_cols = [col for col in df.columns if col not in metadata_cols]

        # 3. Séparation des Cas et des Témoins
        cases = df[df["is_case"] == 1].set_index("match_group_id")[protein_cols]
        controls = df[df["is_case"] == 0]

        # Si k > 1 témoins par cas, on calcule la moyenne des témoins par groupe d'appariement
        controls_mean = controls.groupby("match_group_id")[protein_cols].mean()

        # 4. Calcul vectorisé du Delta P par groupe d'appariement (Cas - Moyenne des Témoins)
        delta_p = cases - controls_mean

        # 5. Calcul des métriques statistiques sur les 2 923 protéines
        mean_series = delta_p.mean(axis=0)  # Moyenne
        var_series = delta_p.var(axis=0)  # Variance
        std_series = delta_p.std(axis=0)  # Écart-type (STD)
        n_samples = len(delta_p)
        sem_series = std_series / np.sqrt(
            n_samples
        )  # Erreur type de la moyenne (SEM)

        # 6. Assemblage du DataFrame de résultats final
        stats_df = pd.DataFrame(
            {
                "mean": mean_series,
                "var": var_series,
                "std": std_series,
                "sem": sem_series,
            }
        )

        stats_df.index.name = "protein_id"
        print(
            f"[RABITPipeline] Calcul terminé pour {len(stats_df)} protéines."
        )

        return stats_df

    def run_full_pipeline(
        self,
        cohort_df: pd.DataFrame,
        work_dir: str,
        output_stats_csv: str,
    ) -> pd.DataFrame:
        """Méthode haut niveau qui enchaîne : export CSV -> inférence GPU ->

        calcul Delta P -> sauvegarde finale.
        """
        # Chemins temporaires pour les échanges de fichiers
        pred_csv = os.path.join(work_dir, "pred_times.csv")
        parquet_dir = os.path.join(work_dir, "rabit_output")

        # 1. Export du fichier CSV pour RABIT
        self.export_prediction_times(cohort_df, pred_csv)

        # 2. Exécution du modèle RABIT sur GPU
        parquet_file = self.run_inference(pred_csv, parquet_dir)

        # 3. Calcul des statistiques protéomiques vectorisées
        stats_df = self.compute_protein_statistics(parquet_file, cohort_df)

        # 4. Sauvegarde du tableau statistique final
        stats_df.to_csv(output_stats_csv)
        print(
            f"[RABITPipeline] Livrable final enregistré sous : {output_stats_csv}"
        )

        return stats_df