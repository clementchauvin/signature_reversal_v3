import femr.datasets
import pandas as pd
from dictionnaire_correspondance_phecode_icd import prepare_icd_to_phecode_mapping
icd_to_phecodes = prepare_icd_to_phecode_mapping()
db = femr.datasets.PatientDatabase("../../../../remote/private/starr_omop_deid/rabit/stanford_all_patients_5main_2025_04_14/extracts/extract") # ou votre chemin STARR
print(len(db))
def list_all_patients_that_have_disease(TARGET_PHECODE = "250.2"):

    cases = []  # Contiendra une liste: {
                #"patient_id": patient_id,
                #"prediction_time": first_diagnosis_date,
                #"gender": patient.gender,
               #"bday": patient_birthday
    i=0
    
    for patient_id in db:
        i=i+1
        patient = db[patient_id]
        first_diagnosis_date = None
        for event in patient.events:
            # Vérifier si l'événement est un diagnostic ICD
            if event.code.startswith("Birth/"):
                            patient_birthday = event.start
                            continue
            if "/" in event.code:
                vocab, code_raw = event.code.split("/", 1)
                code_clean = code_raw.replace(".", "").upper()
            
                # Vérifier si ce code correspond au PheCode cible
                if (vocab, code_clean) in icd_to_phecodes:
                    if any(p.startswith(TARGET_PHECODE) for p in icd_to_phecodes[(vocab, code_clean)]):
                        if first_diagnosis_date is None or event.start < first_diagnosis_date:
                            first_diagnosis_date = event.start
                        
        if first_diagnosis_date is not None:
            cases.append({
                "patient_id": patient_id,
                "prediction_time": first_diagnosis_date,
                "gender": patient.gender,
                "bday": patient_birthday
            })
        if i==1000000:
            break

    cases_df = pd.DataFrame(cases)
    return cases_df

print(len(list_all_patients_that_have_disease()))
