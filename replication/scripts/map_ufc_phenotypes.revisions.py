import numpy as np
import os
import pandas as pd
import pandas_gbq

# Define what individuals we will use in our analysis
v7_samples = np.loadtxt("v7.subjects", dtype=int)
v9_samples = np.loadtxt("v9.subjects", dtype=int)
new_samples = np.loadtxt("new_v8v9.subjects", dtype=int)
ufc_samples = np.loadtxt("ufc_subjects.aou.list", dtype=int)

# Directory where old cohorts are held; we need to see how concordant new cancer phenotypes are to the old ones
old_cohorts = "./"

# Define the 11 cancer types and the corresponding codes we will use 443392: "pancancer", 317510: "leukemia"
omop_codes = {
    4112853: "breast",
    4178976: "thyroid",
    197508: "bladder",
    4038838: "non-hodgkin lymphoma",
    4181351: "ovary",
    197230: "uterus",
    196653: "kidney",
    4162276: "melanoma",
    4163261: "prostate",
    37168850: "colorectal",
    443388: "lung"

}

#cutaneous melanoma, non-medullary thyroid cancer, corpus uteri cancer
omop_codes_exclude = {
	"breast": [4244051,4247348,4003684,4157447,135489],
	"thyroid": [4111011],
	"uterus": [196359],
	"melanoma": [4116085,36712732, 4170619,]

}

def return_query(omop_code: int):
    return f"""
    SELECT
        co.person_id,
        co.condition_concept_id,
        MIN(co.condition_start_date) AS condition_start_date
    FROM `wb-silky-artichoke-2408.C2025Q4R6.condition_occurrence` co
    JOIN `wb-silky-artichoke-2408.C2025Q4R6.concept_ancestor` ca
        ON co.condition_concept_id = ca.descendant_concept_id
    WHERE ca.ancestor_concept_id = {omop_code}
        AND ca.min_levels_of_separation > 0
    GROUP BY co.person_id, co.condition_concept_id
    """

def assign_individuals_to_phenotypes(omop_code: int, revision_dir: str):
    cancer = omop_codes.get(omop_code)
    df = pd.read_gbq(return_query(omop_code))
    df.to_csv(f"{cancer}.tsv", sep="\t", index=False)
    df_ufc = df[df["person_id"].isin(ufc_samples)].copy()
    df_ufc.to_csv(f"{revision_dir}/{cancer}.ufc.tsv", sep="\t", index=False)
    df_v9 = df[df["person_id"].isin(new_samples)].copy()
    df_v9.to_csv(f"{revision_dir}/{cancer}.v9.tsv", sep="\t", index=False)

def compare_aou_phenotypes(cancer: str, old_cohort_dir: str, revision_dir: str):
    ufc_df = pd.read_csv(f"{old_cohort_dir}/{cancer}.metadata", sep="\t")
    ufc_cases = ufc_df.loc[ufc_df["case_control"] == 1, "original_id"].tolist()
    reclassified_ufc_df = pd.read_csv(f"{revision_dir}/{cancer}.ufc.tsv", sep="\t")
    reclassified_cases = reclassified_ufc_df["person_id"].tolist()
    ufc_cases_set = set(ufc_cases)
    reclassified_cases_set = set(reclassified_cases)
    ufc_not_reclassified = ufc_cases_set - reclassified_cases_set
    reclassified_not_ufc = reclassified_cases_set - ufc_cases_set

    print(f"\n{cancer}")
    print("-" * 50)
    print(f"UFC cases:                  {len(ufc_cases_set)}")
    print(f"Reclassified cases:         {len(reclassified_cases_set)}")
    print(f"UFC cases not in reclassified update: {len(ufc_not_reclassified)}")
    print(ufc_not_reclassified)
    # print(f"Reclassified cases not UFC: {len(reclassified_not_ufc)}")
    in_ufc = reclassified_ufc_df["person_id"].isin(ufc_cases_set)

    print("\ncondition_concept_id counts:")
    print("\nCases also in UFC:")
    print(reclassified_ufc_df.loc[in_ufc, "condition_concept_id"].value_counts())
    # print("\nCases NOT in UFC:")
    # print(reclassified_ufc_df.loc[~in_ufc, "condition_concept_id"].value_counts())

    # return {"ufc_not_reclassified": ufc_not_reclassified, "reclassified_not_ufc": reclassified_not_ufc}

# def get_people_with_knowledgeable_about_fhx():
#     query = """
#     SELECT DISTINCT person_id, answer
#     FROM `wb-silky-artichoke-2408.C2025Q4R6.cb_review_survey`
#     WHERE LOWER(question) LIKE "%how much do you know about illnesses or health problems for your parents, grandparents, brothers, sisters, and/or children?%"
#         AND (
#             LOWER(answer) LIKE "%some%"
#             OR LOWER(answer) LIKE "%a lot%"
#             OR LOWER(answer) LIKE "%none at all%"
#         )
#     """
#     df = pd.read_gbq(query)

#     knowledgeable = df[
#         df["answer"].str.contains("Some|A lot", case=False, na=False)
#     ]["person_id"].tolist()

#     none_at_all = df[
#         df["answer"].str.contains("None at all", case=False, na=False)
#     ]["person_id"].tolist()

#     return knowledgeable, none_at_all

def get_people_knowledgeable_about_fhx():
    query = """
    SELECT DISTINCT person_id, answer
    FROM `wb-silky-artichoke-2408.C2025Q4R6.cb_review_survey`
    WHERE LOWER(question) LIKE "%how much do you know about illnesses or health problems for your parents, grandparents, brothers, sisters, and/or children?%"
    """
    df = pd.read_gbq(query)

    knowledgeable = df[
        df["answer"].str.contains("Some|A lot", case=False, na=False)
    ]["person_id"].tolist()

    none_at_all = df[
        df["answer"].str.contains("None at all", case=False, na=False)
    ]["person_id"].tolist()

    answered = df["person_id"].tolist()

    return knowledgeable, none_at_all, answered

def get_people_with_extensive_family_history():
    query = """
    SELECT person_id, answer
    FROM `wb-silky-artichoke-2408.C2025Q4R6.cb_review_survey`
    WHERE LOWER(answer) LIKE "%cancer%"
    	AND LOWER(answer) LIKE "%including%"
	"""

    df = pd.read_gbq(query)
    df = df[~df["answer"].str.contains("Grandparent|Self|non-cancer", na=False)]
    df["family_member"] = df["answer"].str.split("-", n=1).str[1].str.strip()
    family_counts = df.groupby("person_id")["family_member"].nunique().reset_index(name="n_first_degree_relatives")
    extensive_fhx = family_counts[family_counts['n_first_degree_relatives'] >= 3]['person_id'].tolist()
    any_fhx = family_counts[(family_counts['n_first_degree_relatives'] == 1) | (family_counts['n_first_degree_relatives'] == 2)]['person_id'].tolist()

    return extensive_fhx, any_fhx

def compare_fhx_knowledge(extensive_fhx, any_fhx, knowledgeable, none_at_all):
    extensive_fhx = set(extensive_fhx)
    knowledgeable = set(knowledgeable)
    none_at_all = set(none_at_all)

    extensive_knowledgeable = extensive_fhx & knowledgeable
    extensive_none_at_all = extensive_fhx & none_at_all
    neither = extensive_fhx - knowledgeable - none_at_all

    print("\n Extensive Reported Family History of Cancer")
    print(f"Extensive FHx - Count:              {len(extensive_fhx)}")
    print(f"Extensive FHx && Knowledgeable (Some/A lot) - Count: {len(extensive_knowledgeable)}")
    print(f"Extensive FHx && None at all - Count:                {len(extensive_none_at_all)}")
    print(f"Extensive FHx && Did not choose for Knowledgeable or none at all - Count:                    {len(neither)}")
    print(f"Percent of Samples w/ Extensive FHx and are Knowledgeable: {len(extensive_knowledgeable) / len(extensive_fhx) * 100:.1f}%")
    print(f"Percent of Samples w/ Extensive FHx and are Know 'Nothing':   {len(extensive_none_at_all) / len(extensive_fhx) * 100:.1f}%")

    any_fhx = set(any_fhx)
    any_knowledgeable = any_fhx & knowledgeable
    any_none_at_all = any_fhx & none_at_all
    any_neither = any_fhx - knowledgeable - none_at_all

    print("\nPercentage of individuals with Knowledge and no reported FHx")
    print(f"Any FHx - Count: {len(any_fhx)}")
    print(f"Any FHx && Knowledgeable (Some/A lot) - Count: {len(any_knowledgeable)}")
    print(f"Any FHx && None at all - Count: {len(any_none_at_all)}")
    print(f"Any FHx && Neither - Count: {len(any_neither)}")
    print(f"Percent of Samples w/ Any FHx and are Knowledgeable: {len(any_knowledgeable) / len(any_fhx) * 100:.1f}%")
    print(f"Percent of Samples w/ Any FHx and are Know 'Nothing': {len(any_none_at_all) / len(any_fhx) * 100:.1f}%")

    knowledgeable_no_fhx = set(knowledgeable) - set(any_fhx) - set(extensive_fhx)

    print(f"\nKnowledgeable but no reported FHx - Count: {len(knowledgeable_no_fhx)}")
    print(f"Percent of Knowledgeable individuals with no reported FHx: {len(knowledgeable_no_fhx) / len(knowledgeable) * 100:.1f}%")

    answered = set(fhx_knowledge_answered)
    any_fhx = set(any_fhx)
    extensive_fhx = set(extensive_fhx)
    knowledgeable = set(knowledgeable)

    no_reported_fhx = answered - any_fhx - extensive_fhx
    knowledgeable_no_fhx = no_reported_fhx & knowledgeable

    print(f"FHx knowledge question answered: {len(answered)}")
    print(f"No reported FHx: {len(no_reported_fhx)}")
    print(f"Knowledgeable + no reported FHx: {len(knowledgeable_no_fhx)}")
    print(f"Percent of no-reported-FHx who are knowledgeable: {len(knowledgeable_no_fhx) / len(no_reported_fhx) * 100:.1f}%")

if __name__ == "__main__":
    revision_dir = "./revision_cohorts"
    old_cohort_dir = "./analysis"

    extensive_fhx, any_fhx = get_people_with_extensive_family_history()
    #knowledgeable, none_at_all, fhx_knowledge_answered = get_people_knowledgeable_about_fhx()
    #knowledgeable, none_at_all = get_people_with_knowledgeable_about_fhx()

    #compare_fhx_knowledge(extensive_fhx,any_fhx, knowledgeable, none_at_all)

    # for code, cancer in omop_codes.items():
    #     print(f"[INFO]: Assigning Individuals to {cancer}")
    #     assign_individuals_to_phenotypes(code, revision_dir)
    #     print(f"[INFO]: Comparing new vs old {cancer} individuals")
    #     compare_aou_phenotypes(cancer, old_cohort_dir, revision_dir)
    #     break