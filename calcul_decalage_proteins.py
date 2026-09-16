from appariement_malade_sains import match_cohort
import os
import subprocess
import pandas as pd


def decalage_proteins(disease, out_dir="/path/to/your/output_dir"):
    """Calcule les statistiques du décalage protéique (Delta P) entre cas et
    témoins.

    --------------------------------------------------------------------------------
    FORME DU LIVRABLE (Return Value)
    --------------------------------------------------------------------------------
    Objet retourné : pandas.DataFrame (2 923 lignes × 4 colonnes)
    Index          : Nom des protéines (ex: 'INS_prediction', 'TNF_prediction', ...)

    Structure des colonnes :
    +-----------------+-----------+-------------------------------------------------+
    | Nom Colonne     | Type      | Description                                     |
    +-----------------+-----------+-------------------------------------------------+
    | mean_delta      | float64   | Moyenne du décalage (Cas - Témoins) sur N paires|
    | var_delta       | float64   | Variance du décalage entre les paires (s²)      |
    | std_delta       | float64   | Écart-type du décalage (s = sqrt(var_delta))    |
    | sem_delta       | float64   | Erreur-type de la moyenne (s / sqrt(N))         |
    +-----------------+-----------+-------------------------------------------------+

    Aperçu visuel de la sortie :
                        mean_delta  var_delta  std_delta  sem_delta
    INS_prediction        0.420101   0.022500   0.150000   0.003354
    TNF_prediction       -0.120530   0.006400   0.080000   0.001788
    ...                        ...        ...        ...        ...
    --------------------------------------------------------------------------------
    """
    # 1. Obtenir la cohorte appariée
    cohort = match_cohort(disease)

    # 2. Exporter et exécuter RABIT sur GPU
    pred_times_path = os.path.join(out_dir, "t2d_cohort_pred_times.csv")
    os.makedirs(out_dir, exist_ok=True)
    cohort[["patient_id", "prediction_time"]].drop_duplicates().to_csv(
        pred_times_path, index=False
    )

    rabit_cmd = [
        "python",
        "/remote/private/starr_omop_deid/rabit/rabit_pipeline.py",
        "--data_source",
        "shc",
        "--pred_times",
        pred_times_path,
        "--out_dir",
        out_dir,
        "--gpu",
        "0",
    ]
    subprocess.run(rabit_cmd, check=True)

    # 3. Charger les prédictions
    synprot_path = os.path.join(
        out_dir, "t2d_cohort_pred_times_synprot.parquet"
    )
    synprot_df = pd.read_parquet(synprot_path)
    merged_df = cohort.merge(
        synprot_df, left_on="patient_id", right_on="patient_ids"
    )

    protein_cols = [
        col for col in synprot_df.columns if col.endswith("_prediction")
    ]

    # 4. Calcul du Delta P par groupe apparié
    cases_prot = (
        merged_df[merged_df["is_case"] == 1]
        .set_index("match_group_id")[protein_cols]
    )
    controls_mean_prot = (
        merged_df[merged_df["is_case"] == 0]
        .groupby("match_group_id")[protein_cols]
        .mean()
    )

    # DataFrame où chaque ligne est un couple/groupe et chaque colonne une protéine
    delta_proteins_per_group = cases_prot - controls_mean_prot

    # 5. Calcul des statistiques globales sur les 2 923 protéines
    stats_df = pd.DataFrame(
        {
            "mean_delta": delta_proteins_per_group.mean(axis=0),
            "std_delta": delta_proteins_per_group.std(axis=0),
            "var_delta": delta_proteins_per_group.var(axis=0),
            "sem_delta": delta_proteins_per_group.sem(
                axis=0
            ),  # Erreur type de la moyenne
        }
    )

    return stats_df