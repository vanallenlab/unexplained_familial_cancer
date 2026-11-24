import argparse
import pandas as pd
import numpy as np


def filter_cases_controls(df, case_filter, control_filter):
    cases = df.query(case_filter).copy()
    controls = df.query(control_filter).copy() if control_filter else df[df["pancancer_cancer"].fillna(0) == 0].copy()
    return cases, controls


def match_controls_to_cases(
    cancer_type: str,
    df: pd.DataFrame,
    case_filter: str,
    control_filter: str,
    n_matched_controls=30,
    exclude_files=None,
    seed=42,
    output=None
):

    case_age = cancer_type + "_age"

    # Exclude samples included in files designated for removal
    if exclude_files:
        excluded_ids = set()
        for filepath in exclude_files:
            with open(filepath) as f:
                excluded_ids.update(line.strip() for line in f if line.strip())

        df = df[~df["person_id"].astype(str).isin(excluded_ids)]

    # Identify cases and eligible controls
    cases, controls = filter_cases_controls(
        df, case_filter, control_filter
    )

    # Random number used as a tie-breaker
    controls["random_number"] = np.random.default_rng(seed).random(
        len(controls)
    )

    selected_control_ids = []
    available_controls = controls.copy()

    for _, case in cases.iterrows():

        # Match on sex, ancestry, and age (<2 years)
        # First attempt: match within 2 years
        candidate_controls = available_controls[
            (available_controls["sex_at_birth"] == case["sex_at_birth"]) &
            (available_controls["ancestry_pred"] == case["ancestry_pred"]) &
            ((available_controls["age"] - case[case_age]).abs() < 2)
        ].copy()

        # Relax age matching to 5 years if fewer than requested controls are available
        if len(candidate_controls) < n_matched_controls:
            candidate_controls = available_controls[
                (available_controls["sex_at_birth"] == case["sex_at_birth"]) &
                (available_controls["ancestry_pred"] == case["ancestry_pred"]) &
                ((available_controls["age"] - case[case_age]).abs() < 5)
            ].copy()

        candidate_controls["age_difference"] = (
            candidate_controls["age"] - case[case_age]
        ).abs()

        candidate_controls = candidate_controls.sort_values(
            ["age_difference", "random_number"]
        )

        selected = candidate_controls.head(n_matched_controls)

        selected_control_ids.extend(selected["person_id"])

        # Prevent controls from being reused
        available_controls = available_controls.drop(selected.index)

    # Combine cases and selected controls
    selected_ids = set(cases["person_id"]) | set(selected_control_ids)

    df_matched = df[df["person_id"].isin(selected_ids)].copy()

    # Write output
    if output is None:
        output = f"ufc_replication.{cancer_type}.CDRv9.cohort.metadata"

    df_matched.to_csv(output, sep="\t", index=False)

    print(f"Cancer type: {cancer_type}")
    print(f"Cases: {len(cases)}")
    print(f"Matched controls: {len(selected_control_ids)}")
    print(f"Total individuals: {len(df_matched)}")
    print(f"Output: {output}")

    return df_matched


def main():

    parser = argparse.ArgumentParser(description="Match controls to cancer cases.")

    parser.add_argument("-d", "--data", required=True, help="Input metadata TSV")
    parser.add_argument("-c", "--cancer-type", required=True, help="Cancer type (e.g., breast, colon, kidney)")
    parser.add_argument("--control-filter", required=False, default=None, help="Pandas query expression defining eligible controls")
    parser.add_argument("--case-filter", default=None, help="Pandas query expression defining cases. Defaults to '<cancer_type>_cancer == 1'")
    parser.add_argument("-n", "--n-matched-controls", type=int, default=30, help="Maximum number of controls per case (default: 30)")
    parser.add_argument("--exclude-files", nargs="+", default=None, help="Files containing person IDs to exclude")
    parser.add_argument("-o", "--output", default=None, help="Output filename")

    args = parser.parse_args()

    case_filter = args.case_filter
    if case_filter is None:
        case_filter = f"{args.cancer_type}_cancer == 1"

    df = pd.read_csv(args.data, sep="\t", low_memory=False, dtype={"person_id": str})

    match_controls_to_cases(
        cancer_type=args.cancer_type,
        df=df,
        case_filter=case_filter,
        control_filter=args.control_filter,
        n_matched_controls=args.n_matched_controls,
        exclude_files=args.exclude_files,
        seed=42,
        output=args.output
    )


if __name__ == "__main__":
    main()

"""python3 /opt/match_controls.py \
    -d CDRv9.ufc_replication.metadata \
    -c breast \
    --control-filter "pancancer_cancer == 0" \
    -n 30 \
    --seed 42"""

"""python3 /opt/match_controls.py \
    -d CDRv9.ufc_replication.metadata \
    -c breast \
    --control-filter "pancancer_cancer == 0" \
    --exclude-files v7_subjects.list PGV_carriers.list"""