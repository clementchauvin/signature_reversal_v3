
from dictionnaire_correspondance_phecode_icd import prepare_icd_to_phecode_mapping

import femr.datasets
import time
import pandas as pd
exclusion_df = pd.read_csv("original_phecodes_exclusion.csv") 
TARGET_PHECODE = "250.2"
excluded_phecodes = set(
        exclusion_df[exclusion_df['code'] == TARGET_PHECODE]['exclusion_criteria'].astype(str)
    )
# Ajouter la maladie elle-même à l'exclusion
excluded_phecodes.add(TARGET_PHECODE)
# --- 1. À MODIFIER ---
# Chemin direct vers le dossier "extract"
EXTRACT_DIR = "../../../../remote/private/starr_omop_deid/rabit/stanford_all_patients_5main_2025_04_14/extracts/extract"
db = femr.datasets.PatientDatabase(EXTRACT_DIR) # ou votre chemin STARR
icd_to_phecodes = prepare_icd_to_phecode_mapping()
# 1. PRÉ-CALCUL (Exécuté UNE SEULE FOIS en 2 millisecondes)
# On crée un set contenant la chaîne exacte "VOCAB/CODE" issue de FEMR/STARR
EXCLUDED_FEMR_CODES = set()

for (vocab, code_clean), phecodes in icd_to_phecodes.items():
    if not phecodes.isdisjoint(excluded_phecodes):
        # On reconstitue le format exact tel qu'enregistré dans FEMR
        EXCLUDED_FEMR_CODES.add(f"{vocab}/{code_clean}")
        
        # Si FEMR conserve parfois les points dans les codes bruts :
        # On peut ajouter aussi la version avec points si nécessaire

eligible_controls = []  # Contiendra : patient_id, gender
# 2. BOUCLE PATIENT ULTRA-RAPIDE
for patient_id in db:
    
        
    has_exclusion = False
    patient = db[patient_id]
    
    for event in patient.events:
        # Aucun split, aucun replace ! 
        # Une simple vérification d'appartenance dans une table de hachage O(1)
        if event.code in EXCLUDED_FEMR_CODES:
            has_exclusion = True
            break  # Exclu instantanément !
if not has_exclusion:
            eligible_controls.append({
                "patient_id": patient_id,
                "gender": patient.gender
            })

print(f"Nombre de témoins éligibles pour {TARGET_PHECODE}: {len(eligible_controls)}")

