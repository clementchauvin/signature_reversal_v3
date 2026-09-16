import pandas as pd
from collections import defaultdict

#Donne la correspondnace de code ICD10CM ou ICD9CM vers PheCode. Le dictionnaire est de la forme :
# {("ICD10CM", "E119"): {"250.2", "250"}
#Il y a 1876 phecodes uniques. Un code ICD peut correspondre à plusieurs phecodes.
def prepare_icd_to_phecode_mapping(path_to_csv="phecodes_cm.csv"):
    """
    Prépare un dictionnaire de correspondance entre les codes ICD et les PheCodes.
    La clé est un tuple (vocabulaire, code) et la valeur est un ensemble de PheCodes.
    """
# 1. Charger la table de mapping Unrolled
    phecode_df = pd.read_csv(path_to_csv)

    # Nettoyer les codes (supprimer les points et mettre en majuscules)
    phecode_df['clean_icd'] = phecode_df['ICD_CODE'].astype(str).str.replace('.', '', regex=False).str.upper()

    # Dictionnaire : (vocabulaire, code) -> set de PheCodes
    # Exemple de clé : ("ICD10CM", "E119") -> {"250.2", "250"}
    icd_to_phecodes = defaultdict(set)
    for _, row in phecode_df.iterrows():
        key = (row['VOCABULARY_ID'], row['clean_icd'])
        icd_to_phecodes[key].add(str(row['PheCode']))

    return(icd_to_phecodes)

