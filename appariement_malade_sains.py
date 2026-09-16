import random
from collections import defaultdict
import pandas as pd
from has_enough_history import is_eligible_for_history


def match_cohort(
    cases_df,
    controls_df,
    db,
    ponderation_age,
    ponderation_event_count,
    k=4,
    age_tolerance_years=2
   
):
    """Apparie cases_df et controls_df pour générer la cohorte congelée MOTOR."""
    # 1. Dictionnaires d'accès rapide O(1)
    ctrl_bday_dict = dict(zip(controls_df["patient_id"], controls_df["bday"]))
    ctrl_gender_dict = dict(
        zip(controls_df["patient_id"], controls_df["gender"])
    )

    # Indexation des témoins par Sexe
    controls_by_gender = defaultdict(list)
    for p_id, g in ctrl_gender_dict.items():
        controls_by_gender[g].append(p_id)

    available_controls = set(controls_df["patient_id"])
    matched_records = []

    print(
        f"Lancement du matching 1:{k} sur {len(cases_df)} cas identifiés..."
    )

    for _, case in cases_df.iterrows():
        case_id = case["patient_id"]
        t0 = case["prediction_time"]
        case_gender = case["gender"]
        case_bday = case["bday"]

        if pd.isnull(case_bday) or not case_gender:
            continue

        case_age_t0 = (t0 - case_bday).days / 365.25

        # Recherche de candidats témoins
        candidate_ids = controls_by_gender[case_gender]
        scored_candidates = []

        for ctrl_id in candidate_ids:
            if ctrl_id not in available_controls:
                continue

            ctrl_bday = ctrl_bday_dict[ctrl_id]#faut que le témoin ait un birthday pour calculer l'âge à t0
            if pd.isnull(ctrl_bday):
                continue

            # A. Filtrage Âge à t0
            ctrl_age_t0 = (t0 - ctrl_bday).days / 365.25
            age_diff = abs(case_age_t0 - ctrl_age_t0)
            if age_diff > age_tolerance_years:
                continue

            

            if not is_eligible_for_history(db[ctrl_id], db[case_id], t0)[0]:
                continue

            event_diff= is_eligible_for_history(db[ctrl_id], db[case_id], t0)[1]
            score = age_diff * ponderation_age + event_diff * ponderation_event_count
            scored_candidates.append((score, ctrl_id))

        # Sélection des k plus proches voisins
        scored_candidates.sort(key=lambda x: x[0])

        if len(scored_candidates) >= k:
            selected_controls = [
                ctrl_id for _, ctrl_id in scored_candidates[:k]
            ]

            # Cas
            matched_records.append(
                {
                    "patient_id": case_id,
                    "prediction_time": t0,
                    "is_case": 1,
                    "match_group_id": case_id,
                }
            )

            # Témoins appariés (héritent de t0)
            for ctrl_id in selected_controls:
                matched_records.append(
                    {
                        "patient_id": ctrl_id,
                        "prediction_time": t0,
                        "is_case": 0,
                        "match_group_id": case_id,
                    }
                )
                available_controls.remove(ctrl_id)

    final_cohort_df = pd.DataFrame(matched_records)
    print(
        f"Cohorte finalisée : {len(final_cohort_df)} enregistrements prêts."
    )
    return final_cohort_df