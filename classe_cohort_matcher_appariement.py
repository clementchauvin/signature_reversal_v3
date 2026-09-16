import numpy as np
import pandas as pd


class CohortMatcher:
    """Classe responsable de l'appariement 1:k sans remplacement, contrôlant :

    - Le Sexe (Matching exact)
    - L'Âge au moment du diagnostic (+/- max_age_diff_years)
    - La présence d'un suivi médical actif au t0 du cas
    - L'ÉPAISSEUR DU DOSSIER : Nombre de lignes/événements comparables (ratio
    max_event_count_ratio)
    """

    def __init__(
        self,
        k_neighbors: int = 5,
        max_age_diff_years: float = 2.0,
        max_event_count_ratio: float = 2.0,  # Le témoin doit avoir entre (1/ratio) et (ratio) fois le nombre de lignes du cas
    ):
        """Initialise les critères d'appariement.

        :param k_neighbors: Nombre de témoins par cas (ex: 5).
        :param max_age_diff_years: Écart d'âge maximal en années (ex: 2.0 ans).
        :param max_event_count_ratio: Ratio toléré sur le volume de lignes
            (ex: 2.0 signifie entre 50% et 200% des lignes du cas).
        """
        self.k = k_neighbors
        self.max_age_diff_years = max_age_diff_years
        self.max_event_count_ratio = max_event_count_ratio

    def _compute_age_at_index(
        self, birth_datetime: pd.Series, index_date: pd.Series
    ) -> pd.Series:
        """Calcule l'âge exact en années à la date d'index t0."""
        return (
            pd.to_datetime(index_date) - pd.to_datetime(birth_datetime)
        ).dt.days / 365.25

    def match(
        self, cases_df: pd.DataFrame, controls_df: pd.DataFrame
    ) -> pd.DataFrame:
        """Exécute l'appariement 1:k et retourne la cohorte appariée finale.

        :param cases_df: DataFrame issu de CohortExtractor (contenant
            event_count).
        :param controls_df: DataFrame issu de CohortExtractor (contenant
            total_lifetime_events).
        :return: DataFrame apparié prêt pour RABIT (match_group_id, patient_id,
            prediction_time, is_case).
        """
        cases = cases_df.copy()
        controls = controls_df.copy()

        # Standardisation des colonnes au format datetime
        cases["prediction_time"] = pd.to_datetime(cases["prediction_time"])
        cases["birth_datetime"] = pd.to_datetime(cases["birth_datetime"])

        controls["birth_datetime"] = pd.to_datetime(controls["birth_datetime"])
        controls["obs_start"] = pd.to_datetime(
            controls["observation_period_start_date"]
        )
        controls["obs_end"] = pd.to_datetime(
            controls["observation_period_end_date"]
        )

        # Calcul de l'âge des cas au moment de leur t0
        cases["age_at_t0"] = self._compute_age_at_index(
            cases["birth_datetime"], cases["prediction_time"]
        )

        # Ensemble pour verrouiller les témoins déjà sélectionnés (Appariement SANS remplacement)
        used_control_ids = set()
        matched_groups = []

        # =========================================================================
        # BOUCLE D'APPARIEMENT SUR CHAQUE CAS
        # =========================================================================
        for group_id, (_, case_row) in enumerate(cases.iterrows()):

            case_id = case_row["patient_id"]
            t0_date = case_row["prediction_time"]
            gender = case_row["gender"]
            case_age = case_row["age_at_t0"]
            case_events = case_row[
                "event_count"
            ]  # Nombre de lignes du dossier du cas

            # ---------------------------------------------------------------------
            # FILTRE 1 : Exclure les témoins déjà réutilisés
            # ---------------------------------------------------------------------
            available_controls = controls[
                ~controls["patient_id"].isin(used_control_ids)
            ]
            if available_controls.empty:
                break

            # ---------------------------------------------------------------------
            # FILTRE 2 : Sexe identique (Gender matching strict)
            # ---------------------------------------------------------------------
            eligible = available_controls[
                available_controls["gender"] == gender
            ].copy()
            if eligible.empty:
                continue

            # ---------------------------------------------------------------------
            # FILTRE 3 : Présence médicale active au t0 du cas
            # (Le témoin doit être suivi à la date t0 du cas)
            # ---------------------------------------------------------------------
            eligible = eligible[
                (eligible["obs_start"] <= t0_date)
                & (eligible["obs_end"] >= t0_date)
            ]
            if eligible.empty:
                continue

            # ---------------------------------------------------------------------
            # FILTRE 4 : ÉPAISSEUR DU DOSSIER (NOMBRE DE LIGNES / ÉVÉNEMENTS)
            # Le témoin doit avoir un nombre de lignes médicales du même ordre
            # de grandeur que le cas pour éviter d'opposer un dossier riche à un dossier vide
            # ---------------------------------------------------------------------
            min_lines = case_events / self.max_event_count_ratio
            max_lines = case_events * self.max_event_count_ratio

            eligible = eligible[
                (eligible["total_lifetime_events"] >= min_lines)
                & (eligible["total_lifetime_events"] <= max_lines)
            ]
            if eligible.empty:
                continue

            # ---------------------------------------------------------------------
            # FILTRE 5 : Âge similaire au moment du t0 du cas (+/- max_age_diff_years)
            # ---------------------------------------------------------------------
            eligible["age_at_case_t0"] = self._compute_age_at_index(
                eligible["birth_datetime"], t0_date
            )
            eligible["age_diff"] = (
                eligible["age_at_case_t0"] - case_age
            ).abs()

            eligible = eligible[
                eligible["age_diff"] <= self.max_age_diff_years
            ]
            if eligible.empty:
                continue

            # ---------------------------------------------------------------------
            # SÉLECTION : Calcul d'un score de proximité (Combinaison Écart Âge + Écart Volume Lignes)
            # ---------------------------------------------------------------------
            # Différence relative de nombre de lignes
            eligible["lines_diff_ratio"] = (
                (eligible["total_lifetime_events"] - case_events) / case_events
            ).abs()

            # Score combiné (0 = match parfait)
            eligible["match_score"] = (
                eligible["age_diff"] / self.max_age_diff_years
            ) + eligible["lines_diff_ratio"]

            # Trier par le meilleur score et conserver les k premiers témoins
            selected_controls = eligible.sort_values("match_score").head(self.k)

            if len(selected_controls) < 1:
                continue

            # Verrouiller les identifiants sélectionnés
            selected_ids = selected_controls["patient_id"].tolist()
            used_control_ids.update(selected_ids)

            # ---------------------------------------------------------------------
            # CRÉATION DU GROUPE APPARIÉ (1 Cas + k Témoins avec t0 aligné)
            # ---------------------------------------------------------------------
            case_record = pd.DataFrame(
                [
                    {
                        "match_group_id": group_id,
                        "patient_id": case_id,
                        "prediction_time": t0_date,
                        "is_case": 1,
                    }
                ]
            )

            controls_records = pd.DataFrame(
                {
                    "match_group_id": group_id,
                    "patient_id": selected_ids,
                    "prediction_time": t0_date,  # Alignement temporel strict sur le cas !
                    "is_case": 0,
                }
            )

            matched_groups.append(
                pd.concat([case_record, controls_records], ignore_index=True)
            )

        # Assemblage final
        if not matched_groups:
            raise ValueError(
                "L'appariement a échoué : aucun témoin compatible n'a pu être trouvé avec ces critères d'épaisseur de dossier et d'âge."
            )

        return pd.concat(matched_groups, ignore_index=True)