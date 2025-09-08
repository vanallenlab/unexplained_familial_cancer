import yaml
import pandas as pd
from collections import defaultdict

# --- Load YAML file ---
with open('cancer_ontology.yaml', 'r') as f:
    ontology = yaml.safe_load(f)

# Non-Viral Associated Cancers
pgc_cancers = ["appendix","biliary","bladder","bone","brain","breast","colorectal","esophagus","eye","kidney",
"leukemia","lung","lymphoma","melanoma","meninges","myelomastocytic","nervous","neuroendocrine","ovary","pancreas","parathyroid",
"prostate","stomach","small intestines","testis","thyroid","uterus"]

cancel_out = {"intrahepatic_bile_duct_carcinoma":"Liver","primary_malignant_neoplasm_of_renal_pelvis":"Kidney",
"neuroendocrine_carcinoma_of_appendix":"Appendix","malignant_neuroendocrine_tumor_of_duodenum":"Small_Intestines",
"malignant_neuroendocrine_tumor_of_ileum":"Small_Intestines", "malignant_neuroendocrine_tumor_of_small_intestine":"Small_Intestines",
"primary_malignant_neuroendocrine_neoplasm_of_duodenum":"Small_Intestines","primary_malignant_neuroendocrine_neoplasm_of_ileum":"Small_Intestines",
"primary_malignant_neuroendocrine_neoplasm_of_small_intestine":"Small_Intestines","primary_malignant_neoplasm_of_islets_of_langerhans":"Pancreas",
"malignant_carcinoid_tumor_of_lung":"Lung","malignant_carcinoid_tumor_of_bronchus":"Lung","malignant_carcinoid_tumor_of_small_intestine":"Small_Intestines",
"malignant_carcinoid_tumor_of_ileum":"Small_Intestines","malignant_carcinoid_tumor_of_ileum":"Small_Intestines","malignant_carcinoid_tumor_of_colon":"Colorectal",
"carcinoid_tumor_of_small_intestine":"Small_Intestines","carcinoid_tumor_of_large_intestine":"Colorectal","carcinoid_tumor_of_intestine":"Small_Intestines",
"carcinoid_tumor_of_ileum":"Small_Intestines","carcinoid_tumor_of_stomach":"Stomach","carcinoid_tumor_of_lung":"Lung","burkitts_lymphoma_clinical":"Non-Hodgkins",
"primary_malignant_neoplasm_of_urethra":"Prostate","mesothelioma_malignant,_clinical_disorder":"Lung","malignant_mesothelioma_of_pleura":"Lung",
"primary_malignant_neoplasm_of_pleura":"Lung","primary_malignant_neoplasm_of_pleura":"Lung","large_cell_anaplastic_lymphoma":"Non-Hodgkins",
"extranodal_marginal_zone_bcell_lymphoma_of_mucosaassociated_lymphoid_tissue_maltlymphoma":"Non-Hodgkins"}

# --- Build a flat mapping from diagnosis to system labels ---
dx_to_systems = {}

systems = ontology['cancer_ontology']['systems']
for system in systems:
    sys_name = system['name'].replace(' ', '_')

    # Top-level types
    for dx in system.get('types', []):
        dx_to_systems[dx] = [sys_name]
    # Sub-systems
    for sub in system.get('sub_systems', []):
        sub_name = sub['name'].replace(' ', '_')
        for dx in sub.get('types', []):
            dx_to_systems[dx] = [sub_name, sys_name]

# --- Load TSV file ---
df = pd.read_csv('dfci-ufc.aou.phenos.tsv', sep='\t')

# --- Grouping function ---
def group_dxs(row):
    if row.get('cancer', '').strip().lower() == 'control':
        return 'control'

    original_dx = str(row.get('original_dx', '')).strip()
    if not original_dx or original_dx.lower() == 'nan':
        return 'NA'

    systems_found = set()

    for dx in original_dx.split(';'):
        dx = dx.strip()
        if dx in dx_to_systems:
            systems_found.update(dx_to_systems[dx])
            #systems_found.add(';'.join(dx_to_systems[dx]))

    if systems_found:
        return ';'.join(sorted(systems_found))
    else:
        return 'NA'

# --- Apply function to create new column ---
df['original_dx_grouped'] = df.apply(group_dxs, axis=1)

def clean_dx_grouped(row):
    labels = set(row['original_dx_grouped'].split(';'))

    # Remove any sub-systems based on cancel_out keys
    for key, value in cancel_out.items():
        if key in str(row['original_dx']).lower() and value in labels:
            labels.discard(value)

    # Build parent → children mapping from dx_to_systems
    parent_to_children = {}
    for systems in dx_to_systems.values():
        if len(systems) > 1:
            child, parent = systems[0], systems[1]
            parent_to_children.setdefault(parent, set()).add(child)

    # Remove parent if none of its children remain
    for parent, children in parent_to_children.items():
        # intersect with current labels to see if any child exists
        if parent in labels and len(labels.intersection(children)) == 0:
            labels.discard(parent)

    return ';'.join(sorted(labels))

# Apply to dataframe
df['original_dx_grouped'] = df.apply(clean_dx_grouped, axis=1)

# --- Function to count PGC cancers per patient (assume primary unless ambiguous system) ---
def count_pgc(row):
    if row.get('cancer', '').strip().lower() == 'control':
        return 0

    original_dx = str(row.get('original_dx', '')).strip().lower()
    if not original_dx or original_dx == 'nan':
        return 0

    ambiguous_systems = {"lung", "bone", "brain", "liver"}

    # Split diagnoses by semicolon
    diagnoses = [dx.strip() for dx in original_dx.split(';')]

    mapped_systems = []
    for dx in diagnoses:
        if dx in dx_to_systems:
            systems = dx_to_systems[dx]
            for sys in systems:
                sys_lower = sys.lower()
                # Keep if non-ambiguous OR explicitly marked as "primary"
                if sys_lower not in ambiguous_systems or "primary" in dx:
                    mapped_systems.append(sys)

    # Deduplicate to avoid double counting
    unique_systems = set(sys.lower() for sys in mapped_systems)

    # Count overlap with PGC cancers
    return sum(sys in [c.lower() for c in pgc_cancers] for sys in unique_systems)

# Apply to dataframe
df['Possibly_Genetic_Cancers'] = df.apply(count_pgc, axis=1)

# --- Save result ---
df.to_csv('dfci-ufc.aou.phenos.v2.tsv', sep='\t', index=False)
