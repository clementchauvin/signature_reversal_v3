def get_patient_history_stats(patient, t0):
    """
    Calcule la profondeur d'historique (en jours) et le nombre d'événements avant t0.
    """
    events_before_t0 = [e for e in patient.events if e.start <= t0]
    
    if not events_before_t0:
        return 0, 0  # Aucun historique avant t0
        
    first_event_date = min(e.start for e in events_before_t0)
    history_depth_days = (t0 - first_event_date).days
    event_count = len(events_before_t0)
    
    return history_depth_days, event_count