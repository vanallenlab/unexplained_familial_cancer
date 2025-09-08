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
    "extranodal_marginal_zone_bcell_lymphoma_of_mucosaassociated_lymphoid_tissue_maltlymphoma": "Lymphoma"
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
    "appendix","biliary","bladder","bone","brain","breast","colorectal","esophagus","eye","kidney",
    "leukemia","lung","lymphoma","melanoma","meninges","myelomastocytic","nervous","neuroendocrine",
    "ovary","pancreas","parathyroid","prostate","stomach","small_intestines","soft_tissue",
    "testis","thyroid","uterus"
]

# --- Count function ---
def count_pgc(row):
    if row.get("cancer", "").strip().lower() == "control":
        return 0

    grouped = str(row.get("original_dx_grouped", "")).strip()
    if not grouped or grouped.lower() == "nan":
        return 0

    ambiguous = {"lung","liver","bone","brain","meninges"}
    dxs = str(row.get("original_dx", "")).lower()

    systems_found = set(grouped.split(";"))
    filtered = set()

    for sys in systems_found:
        sys_lower = sys.lower()
        if sys_lower in ambiguous:
            if "primary" in dxs:
                filtered.add(sys_lower)
        else:
            filtered.add(sys_lower)

    count = sum(sys in [c.lower() for c in pgc_cancers] for sys in filtered)

    # bump ambiguous-only cases
    if count == 0:
        for sys in systems_found:
            if sys.lower() in ambiguous:
                return 1
    return count

# --- Main ---
if __name__ == "__main__":
    df = pd.read_csv("dfci-ufc.aou.phenos.tsv", sep="\t")
    df["original_dx_grouped"] = df.apply(group_dxs, axis=1)
    df["Possibly_Genetic_Cancers"] = df.apply(count_pgc, axis=1)
    df.to_csv("dfci-ufc.aou.phenos.v2.tsv", sep="\t", index=False)

