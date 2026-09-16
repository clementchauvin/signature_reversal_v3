import femr.datasets


# The database contains all the features used for the patients
database = femr.datasets.PatientDatabase('../../../../remote/private/starr_omop_deid/rabit/stanford_all_patients_5main_2025_04_14/extracts/extract')
#database = femr.datasets.PatientDatabase('../../synthetic_data/extract_lite')

patients = list(database)


# For example, we can load a single patient with patient_id 110
for i in range(20000):
    

    patient = database[patients[i]]

    # Print out all the events for the patients
    for event in patient.events:
        if 'ICD' in event.code:# and event.value is not None:
            print(event)