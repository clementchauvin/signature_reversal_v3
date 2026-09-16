from datetime import datetime, timedelta
from get_patient_history_length import get_patient_history_stats

def is_eligible_for_history(patient_sain, patient_malade, t0, tolerance_ratio=0.25):
    """
    Vérifie si le patient a une profondeur d'historique suffisante avant t0.
    """
    depth_days_sain, event_count_sain = get_patient_history_stats(patient_sain, t0)
    depth_days_malade, event_count_malade = get_patient_history_stats(patient_malade, t0)

   
    return((event_count_sain >= (1-tolerance_ratio) * event_count_malade and event_count_sain <= event_count_malade * (1+tolerance_ratio)), abs(event_count_malade - event_count_sain) <= tolerance_ratio * event_count_malade)

