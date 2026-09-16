import femr.datasets
from collections import defaultdict
from dictionnaire_correspondance_phecode_icd import prepare_icd_to_phecode_mapping


# 1. Charger la base
database = femr.datasets.PatientDatabase('../../synthetic_data/extract_lite') # ou votre chemin STARR


def dictionnaire_patients_maladie(icd_to_phecodes):#icd_to_phecodes is the dictionary from dictionnaire_correspondance_phecode_icd.py
    """
    Crée un dictionnaire qui mappe chaque patient à ses maladies (PheCodes) avec la date de première apparition.
    La structure est la suivante :
    patient_history[pid][phecode] = date_premiere_apparition
    """

    # Structures de stockage en mémoire
    # patient_history[pid][phecode] = date_premiere_apparition
    patient_history = defaultdict(dict)
    patient_birthdays = {}

    for pid in database:
        patient = database[pid]
        for event in patient.events:
            # Récupérer la date de naissance pour l'appariement futur
            if event.code.startswith("Birth/"):
                patient_birthdays[pid] = event.start
                continue
            
            if "/" not in event.code:
                continue
            
            vocab, raw_code = event.code.split("/", 1)
            clean_code = raw_code.replace(".", "").upper()
            key = (vocab, clean_code)
        
            # Si le code ICD est dans le mapping PheCode
            if key in icd_to_phecodes:
                for pcode in icd_to_phecodes[key]:
                    # On ne garde que la date la plus ancienne (T0)
                    if pcode not in patient_history[pid] or event.start < patient_history[pid][pcode]:
                        patient_history[pid][pcode] = event.start

    print(f"Historique indexé pour {len(patient_history)} patients.")
    return patient_history, patient_birthdays
icd_to_phecodes = prepare_icd_to_phecode_mapping()
dictionnaire_patients_maladie(icd_to_phecodes)
  
