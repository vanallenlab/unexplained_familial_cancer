import pandas as pd
import numpy as np

def filter_cases_controls(df, case_filter, control_filter):
    cases = df.query(case_filter).copy()
    controls = df.query(control_filter).copy()
    return cases, controls

def match_cases_to_controls(cancer_type: str, df: pd.DataFrame, case_filter: str, control_filter:str , n_matched_controls=30, exclude_files=None, seed=42):

    case_age = cancer_type + "_age"

    # Exclude samples included in a file designated for removal (e.g., indiviuduals with a PGV, or in CDRv7)
    if exclude_files:
        excluded_ids = {
            line.strip()
            for filepath in exclude_files
            for line in open(filepath)
            if line.strip()
        }
        df = df[~df["person_id"].astype(str).isin(excluded_ids)].copy()

    # Use the filter_cases_controls() function to obtain a suitable set of cases and controls
    cases, controls = filter_cases_controls(df, case_filter, control_filter)

    # Random number will be a tiebreaker for matching controls
    controls["random_number"] = np.random.default_rng(seed).random(len(controls))

    selected_control_ids = []
    available_controls = controls.copy()

    for _, case in cases.iterrows():

        # The only strict requirement for control matching is sex_at_birth and ancestry_pred
        candidate_controls = available_controls[
            (available_controls["sex_at_birth"] == case["sex_at_birth"]) &
            (available_controls["ancestry_pred"] == case["ancestry_pred"]) &
            ((available_controls["age"] - case[case_age]).abs() < 2)
        ].copy()

        selected = candidate_controls.head(n_matched_controls).copy()
        
        selected_control_ids.extend(selected["person_id"])
        available_controls = available_controls.drop(selected.index)

    selected_ids = set(cases["person_id"]) | set(selected_control_ids)

    return df[df["person_id"].isin(selected_ids)].copy()