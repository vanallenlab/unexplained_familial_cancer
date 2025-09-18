import yaml
import pandas as pd
from collections import defaultdict

# --- Load ontology ---
with open("cancer_ontology.v5.yaml", "r") as f:
    ontology = yaml.safe_load(f)

# --- Cancel out dictionary (exact matches only) ---
cancel_out = {
    "intrahepatic_bile_duct_carcinoma": "Liver",
    "liver_cell_carcinoma":"Gastrointestinal",
    "malignant_neoplasm_of_liver":"Gastrointestinal",
    "neoplasm_of_liver":"Gastrointestinal",
    "primary_malignant_neoplasm_of_liver":"Gastrointestinal",
    "primary_malignant_neoplasm_of_renal_pelvis": "Kidney",
    "neuroendocrine_carcinoma_of_appendix": "Appendix",
    "carcinoid_tumor_of_appendix": "Appendix",
    "malignant_neuroendocrine_tumor_of_duodenum": "Small_Intestines",
    "malignant_neuroendocrine_tumor_of_ileum": "Small_Intestines",
    "malignant_neuroendocrine_tumor_of_small_intestine": "Small_Intestines",
    "primary_malignant_neuroendocrine_neoplasm_of_duodenum": "Small_Intestines",
    "primary_malignant_neuroendocrine_neoplasm_of_ileum": "Small_Intestines",
    "primary_malignant_neuroendocrine_neoplasm_of_small_intestine": "Small_Intestines",
    "primary_malignant_neoplasm_of_islets_of_langerhans": "Pancreas",
    "malignant_carcinoid_tumor_of_lung": "Lung",
    "malignant_carcinoid_tumor_of_bronchus": "Lung",
    "malignant_carcinoid_tumor_of_small_intestine": "Small_Intestines",
    "malignant_carcinoid_tumor_of_ileum": "Small_Intestines",
    "malignant_carcinoid_tumor_of_colon": "Colorectal",
    "carcinoid_tumor_of_small_intestine": "Small_Intestines",
    "carcinoid_tumor_of_large_intestine": "Colorectal",
    "carcinoid_tumor_of_intestine": "Small_Intestines",
    "carcinoid_tumor_of_ileum": "Small_Intestines",
    "carcinoid_tumor_of_stomach": "Stomach",
    "carcinoid_tumor_of_lung": "Lung",
    "burkitts_lymphoma_clinical": "Non-Hodgkins",
    "burkitts_lymphoma_clinical": "Lymphoma",
    "primary_malignant_neoplasm_of_urethra": "Prostate",
    "mesothelioma_malignant,_clinical_disorder": "Lung",
    "malignant_mesothelioma_of_pleura": "Lung",
    "primary_malignant_neoplasm_of_pleura": "Lung",
    "large_cell_anaplastic_lymphoma": "Non-Hodgkins",
    "large_cell_anaplastic_lymphoma": "Lymphoma",
    "extranodal_marginal_zone_bcell_lymphoma_of_mucosaassociated_lymphoid_tissue_maltlymphoma": "Non-Hodgkins",
    "extranodal_marginal_zone_bcell_lymphoma_of_mucosaassociated_lymphoid_tissue_maltlymphoma": "Lymphoma",
    "malignant_melanoma_of_right_choroid":"Eye",
    "malignant_melanoma_of_right_choroid":"Nervous"
}

# --- Build dx → systems map + parent tracking ---
dx_to_systems = {}
parent_map = defaultdict(set)  # parent -> set(children)

def add_mapping(dx, systems):
    dx_to_systems[dx] = systems
    if len(systems) > 1:
        child, parent = systems[0], systems[-1]
        parent_map[parent].add(child)

systems = ontology["cancer_ontology"]["systems"]
for system in systems:
    sys_name = system["name"].replace(" ", "_")

    for dx in system.get("types", []):
        add_mapping(dx, [sys_name])

    for sub in system.get("sub_systems", []):
        sub_name = sub["name"].replace(" ", "_")
        for dx in sub.get("types", []):
            add_mapping(dx, [sub_name, sys_name])

def group_dxs(row):
    if row.get('cancer', '').strip().lower() == 'control':
        return 'control'

    original_dx = str(row.get('original_dx', '')).strip()
    if not original_dx or original_dx.lower() == 'nan':
        return 'NA'

    systems_found = set()
    dx_to_keep = []

    for dx in original_dx.split(';'):
        dx = dx.strip()
        if dx in dx_to_systems:
            mapped = dx_to_systems[dx]
            dx_to_keep.append((dx, mapped))
            systems_found.update(mapped)

    # Apply cancel_out
    for dx, mapped in dx_to_keep:
        if dx in cancel_out:
            to_remove = cancel_out[dx]
            if to_remove in systems_found:
                systems_found.remove(to_remove)

    if systems_found:
        return ';'.join(sorted(systems_found))
    else:
        return 'NA'


# --- PGC cancers ---
pgc_cancers = [
    "appendix","anet","biliary","bladder","bone","brain","breast","colorectal","esophagus","eye","gnet","kidney",
    "leukemia","lung","lnet","lymphoma","melanoma","meninges","myelomastocytic","nervous","neuroendocrine",
    "ovary","pancreas","parathyroid","pituitary","pnet","pheochromocytoma","prostate","stomach","small_intestines","soft_tissue",
    "testis","thyroid","thymus","uterus"
]

def count_pgc_direct(row):
    """
    Returns the number of Possibly Genetic Cancers (PGC) for a patient.
    Counts only diagnoses in pgc_cancers.
    Ambiguous cancers (lung, brain, meninges, bone, liver) are only counted if 'primary' appears in the dx.
    """
    if row.get('cancer', '').strip().lower() == 'control':
        return 0

    original_dx = str(row.get('original_dx', '')).strip()
    if not original_dx or original_dx.lower() == 'nan':
        return 0

    # Systems already grouped
    grouped_dx = str(row.get('original_dx_grouped', '')).strip()
    if not grouped_dx or grouped_dx.lower() == 'nan':
        return 0

    # Ambiguous systems that require 'primary' to count
    ambiguous = {"lung","brain","meninges","bone","liver"}

    # Count how many PGC systems are in the grouped dx
    systems_found = set(grouped_dx.split(';'))
    count = 0

    for sys in systems_found:
        sys_lower = sys.lower()
        if sys_lower not in pgc_cancers:
            continue

        if sys_lower in ambiguous:
            # Check if any diagnosis referring to this system contains 'primary'
            for dx in original_dx.split(';'):
                dx = dx.strip().lower()
                if sys_lower in dx and 'primary' in dx:
                    count += 1
                    break
        else:
            count += 1

    # Parent-child additions
    parent_children = {
        "Gastrointestinal": {"Neuroendocrine","Liver","Colorectal","Small_Intestines","Stomach","Appendix","Biliary","Esophagus","Pancreas"},
        "Endocrine": {"Thyroid","Parathyroid"},
        "Urinary": {"Kidney","Bladder"},
        "Nervous": {"Eye","Brain","Meninges"}
    }

    for parent, children in parent_children.items():
        if parent in systems_found and not systems_found.intersection(children):
            count += 1

    return count



def map_and_concat(row):
    terms = set()
    # --- mapping dictionary ---
    mapping = {
        "Bladder": "Bladder",
        "Hemotologic_System": "Blood_Soft_Tissue",
        "Sarcoma": "Blood_Soft_Tissue",
        "Lymphoma": "Blood_Soft_Tissue",
        "Bone": "Bone",
        "Brain": "Brain",
        "Meninges": "Brain",
        "Breast": "Breast",
        "Cervix": "Cervix",
        "Colorectal": "Colorectal",
        "Esophagus": "Esophagus",
        "Eye": "Eye",
        "Kidney": "Kidney",
        "Lung": "Lung",
        "Ovary": "Ovary",
        "Pancreas": "Pancreas",
        "Prostate": "Prostate",
        "Stomach": "Stomach",
        "Thyroid": "Thyroid",
        "Uterus": "Uterus",
        "Oropharynx": "HNSCC",
        "ENT": "HNSCC"
    }
    # --- handle original_dx ---
    if pd.notna(row["original_dx"]):
        for term in str(row["original_dx"]).split(";"):
            term = term.strip()
            if term in mapping:
                terms.add(mapping[term])
            else:
                terms.add(term)  # keep unmapped terms too

    # --- handle cancers_in_FDRs ---
    if pd.notna(row["cancers_in_FDRs"]):
        for term in str(row["cancers_in_FDRs"]).split(";"):
            term = term.strip()
            if term in mapping:
                terms.add(mapping[term])
            else:
                terms.add(term)

    # --- return unique semi-colon joined string ---
    return ";".join(sorted(terms))


# --- Main ---
if __name__ == "__main__":
    df = pd.read_csv("dfci-ufc.aou.phenos.tsv", sep="\t")
    df2 = pd.read_csv("dfci-ufc.aou.family_cancers.tsv", sep='\t',comment="#")
    df["original_dx_grouped"] = df.apply(group_dxs, axis=1)
    df["Possibly_Syndromic_Cancers"] = df.apply(count_pgc_direct, axis=1)
    df = df.merge(df2, left_on = "Sample", right_on = "sample_id",how="left")
    df["all_fam_dx"] = df.apply(map_and_concat, axis=1)
    df = df.fillna("NA")
    df.to_csv("dfci-ufc.aou.phenos.v2.tsv", sep="\t", index=False)

