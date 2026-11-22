
import statsmodels.api as sm
import pandas as pd
import numpy as np


def run_cancer_carrier_logistic(
    carrier_ids_file: str,
    gene_tested: str,
    cancer_type: str,
    dataset: str,
    eur_only: bool = False,
    haplotype_specific_ids: str = None,
    adjust_for_sex: bool = True,
    adjust_for_ancestry: bool = True,
    pgv_carriers: str = None,
    controls_noFHx: bool = True,
    cases_FHx: int = 0,
    pgs_code: str = None,
    provided_output_file: str = None
):

    cancer_col = f"{cancer_type}_cancer"
    output_file = provided_output_file or f"{gene_tested}_{cancer_type}_logreg.stats"

    # Load carrier IDs
    carrier_ids = set(pd.read_csv(carrier_ids_file, header=None, dtype=str)[0].str.strip())

    # Load dataset
    df = pd.read_csv(dataset, sep="\t", index_col=False, dtype={"person_id": str})
    df["person_id"] = df["person_id"].str.strip()

    # Restrict to haplotype-specific individuals if provided
    if haplotype_specific_ids is not None:
        haplotype_ids = set(pd.read_csv(haplotype_specific_ids, header=None, dtype=str)[0].str.strip())
        df = df[df["person_id"].isin(haplotype_ids)]

    # Exclude individuals with PGVs if provided
    if pgv_carriers is not None:
        pgv_carrier_ids = set(pd.read_csv(pgv_carriers, header=None, dtype=str)[0].str.strip())
        df = df[~df["person_id"].isin(pgv_carrier_ids)]

    # Restrict controls to individuals with no reported family history
    if controls_noFHx:
        df = df[
            (df["pancancer_cancer"].fillna(0) == 1) |
            ((df["pancancer_cancer"].fillna(0) == 0) & (df["num_relatives_w_cancer"].fillna(0) == 0))
        ]

    # Restrict cases to the specified minimum family history burden
    df = df[
        ((df["pancancer_cancer"].fillna(0) == 1) & (df["num_relatives_w_cancer"].fillna(0) >= cases_FHx)) |
        (df["pancancer_cancer"].fillna(0) == 0)
    ]

    # Select ancestry and corresponding PCs
    if eur_only:
        df = df[df["ancestry_pred"] == "eur"]
        pcs = [f"EUR_PC{i}" for i in range(1, 6)]
    else:
        pcs = [f"PC{i}" for i in range(1, 6)]

    # Restrict sex for sex-specific cancers
    if cancer_type in ["breast", "ovarian", "endometrial"]:
        df = df[df["sex_at_birth"] == "Female"]
    elif cancer_type in ["prostate", "male_breast"]:
        df = df[df["sex_at_birth"] == "Male"]

    # Define cases and carrier status
    df["case"] = (df[cancer_col] == 1).astype(int)
    df["carrier"] = df["person_id"].isin(carrier_ids).astype(int)

    # Build covariate list
    covariates = ["carrier"]

    if cancer_type not in ["breast", "ovarian", "endometrial", "prostate", "male_breast"] and adjust_for_sex:
        df["sex_binary"] = df["sex_at_birth"].map({"Female": 0, "Male": 1})
        covariates.append("sex_binary")

    if adjust_for_ancestry:
        covariates += pcs

    if pgs_code is not None:
        if pgs_code not in df.columns:
            raise ValueError(f"PGS column '{pgs_code}' not found in dataframe.")
        covariates.append(pgs_code)

    # Remove missing and infinite model values
    model_df = df[["case"] + covariates].replace([np.inf, -np.inf], np.nan).dropna()

    # Fit logistic regression
    X = sm.add_constant(model_df[covariates])
    y = model_df["case"]
    model = sm.Logit(y, X).fit(disp=False)

    # Collect results
    result = {
        "gene_tested": gene_tested,
        "cancer_type": cancer_type,
        "test": "log-reg",
        "num_carriers": int(model_df["carrier"].sum()),
        "num_noncarriers": int((model_df["carrier"] == 0).sum()),
        "eur_only": eur_only,
        "haplotype_provided": haplotype_specific_ids is not None,
        "pgv_excluded": pgv_carriers is not None,
        "controls_noFHx": controls_noFHx,
        "cases_FHx": cases_FHx,
        "sex_adjustment": adjust_for_sex,
        "ancestry_adjustment": adjust_for_ancestry,
        "pgs_code": pgs_code,
        "cases_tested": int(y.sum()),
        "controls_tested": int((y == 0).sum()),
        "OR": np.exp(model.params["carrier"]),
        "P": model.pvalues["carrier"]
    }

    # Save results as a new TSV
    pd.DataFrame([result]).to_csv(output_file, sep="\t", index=False)

    return


def run_cancer_carrier_fisher(
    carrier_ids_file: str,
    gene_tested: str,
    cancer_type: str,
    dataset: str,
    eur_only: bool = False,
    haplotype_specific_ids: str = None,
    pgv_carriers: str = None,
    controls_noFHx: bool = True,
    cases_FHx: int = 0,
    provided_output_file: str = None
):

    cancer_col = f"{cancer_type}_cancer"
    output_file = provided_output_file or f"{gene_tested}_{cancer_type}_fisher.stats"

    # Load carrier IDs
    carrier_ids = set(pd.read_csv(carrier_ids_file, header=None, dtype=str)[0].str.strip())

    # Load dataset
    df = pd.read_csv(dataset, sep="\t", index_col=False, dtype={"person_id": str})
    df["person_id"] = df["person_id"].str.strip()

    # Restrict to haplotype-specific individuals if provided
    if haplotype_specific_ids is not None:
        haplotype_ids = set(pd.read_csv(haplotype_specific_ids, header=None, dtype=str)[0].str.strip())
        df = df[df["person_id"].isin(haplotype_ids)]

    # Exclude individuals with PGVs if provided
    if pgv_carriers is not None:
        pgv_carrier_ids = set(pd.read_csv(pgv_carriers, header=None, dtype=str)[0].str.strip())
        df = df[~df["person_id"].isin(pgv_carrier_ids)]

    # Restrict controls to individuals with no reported family history
    if controls_noFHx:
        df = df[
            (df["pancancer_cancer"].fillna(0) == 1) |
            ((df["pancancer_cancer"].fillna(0) == 0) & (df["num_relatives_w_cancer"].fillna(0) == 0))
        ]

    # Restrict cases to the specified minimum family history burden
    df = df[
        ((df["pancancer_cancer"].fillna(0) == 1) & (df["num_relatives_w_cancer"].fillna(0) >= cases_FHx)) |
        (df["pancancer_cancer"].fillna(0) == 0)
    ]

    # Select ancestry
    if eur_only:
        df = df[df["ancestry_pred"] == "eur"]

    # Restrict sex for sex-specific cancers
    if cancer_type in ["breast", "ovarian", "endometrial"]:
        df = df[df["sex_at_birth"] == "Female"]
    elif cancer_type in ["prostate", "male_breast"]:
        df = df[df["sex_at_birth"] == "Male"]

    # Define cases and carrier status
    df["case"] = (df[cancer_col] == 1).astype(int)
    df["carrier"] = df["person_id"].isin(carrier_ids).astype(int)

    # Construct 2x2 contingency table
    cases_carriers = ((df["case"] == 1) & (df["carrier"] == 1)).sum()
    cases_noncarriers = ((df["case"] == 1) & (df["carrier"] == 0)).sum()
    controls_carriers = ((df["case"] == 0) & (df["carrier"] == 1)).sum()
    controls_noncarriers = ((df["case"] == 0) & (df["carrier"] == 0)).sum()

    contingency_table = [
        [cases_carriers, cases_noncarriers],
        [controls_carriers, controls_noncarriers]
    ]

    OR, P = fisher_exact(contingency_table, alternative="two-sided")

    # Collect results
    result = {
        "gene_tested": gene_tested,
        "cancer_type": cancer_type,
        "test": "fisher",
        "num_carriers": int(df["carrier"].sum()),
        "num_noncarriers": int((df["carrier"] == 0).sum()),
        "eur_only": eur_only,
        "haplotype_provided": haplotype_specific_ids is not None,
        "pgv_excluded": pgv_carriers is not None,
        "controls_noFHx": controls_noFHx,
        "cases_FHx": cases_FHx,
        "cases_tested": int(df["case"].sum()),
        "controls_tested": int((df["case"] == 0).sum()),
        "case_carriers": int(cases_carriers),
        "control_carriers": int(controls_carriers),
        "OR": OR,
        "P": P
    }

    # Save results as a new TSV
    pd.DataFrame([result]).to_csv(output_file, sep="\t", index=False)