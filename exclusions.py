# 1. Charger la table d'exclusion officielle
import pandas as pd
import femr.datasets
exclusion_df = pd.read_csv("original_phecodes_exclusion.csv") 
db = femr.datasets.PatientDatabase("../../../../remote/private/starr_omop_deid/rabit/stanford_all_patients_5main_2025_04_14/extracts/extract") # ou votre chemin STARR


def list_eligible_controls(TARGET_PHECODE, cases_df, icd_to_phecodes):
    # Récupérer tous les PheCodes exclus pour TARGET_PHECODE
    excluded_phecodes = set(
        exclusion_df[exclusion_df['code'] == TARGET_PHECODE]['exclusion_criteria'].astype(str)
    )
    # Ajouter la maladie elle-même à l'exclusion
    excluded_phecodes.add(TARGET_PHECODE)

    # 2. Filtrer les témoins éligibles
    case_ids = set(cases_df['patient_id'])
    eligible_controls = []

    for patient_id in db:
        if patient_id in case_ids:
            continue  # Ce patient est déjà un cas
            
        has_exclusion = False
        patient = db[patient_id]
        
        for event in patient.events:
            if event.code.startswith("Birth/"):
                patient_birthday = event.start
                continue
            if "/" in event.code:
                vocab, code_raw = event.code.split("/", 1)
                code_clean = code_raw.replace(".", "").upper()
                
                matched_phecodes = icd_to_phecodes.get((vocab, code_clean), set())
                if not matched_phecodes.isdisjoint(excluded_phecodes):
                    has_exclusion = True
                    break  # Écarté du groupe témoin
                    
        if not has_exclusion:
            eligible_controls.append({
                "patient_id": patient_id,
                "gender": patient.gender,
                "bday": patient_birthday
            })

        eligible_controls_df = pd.DataFrame(eligible_controls)

    return eligible_controls_df

