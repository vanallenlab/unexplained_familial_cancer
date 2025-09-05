import yaml
import pandas as pd

# --- Load YAML file ---
with open('cancer_ontology.yaml', 'r') as f:
    ontology = yaml.safe_load(f)

# Non-Viral Associated Cancers
pgc_cancers = ["appendix","biliary","bladder","bone","brain","breast","colorectal","esophagus","eye","kidney",
"leukemia","lung","lymphoma","melanoma","meninges","myelomastocytic","nervous","neuroendocrine","ovary","pancreas","parathyroid",
"prostate","stomach","small intestines","testis","thyroid","uterus"]

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

# --- Function to count PGC cancers per patient (based on 'primary' diagnoses in original_dx) ---
def count_pgc(row):
    if row.get('cancer', '').strip().lower() == 'control':
        return 0

    original_dx = str(row.get('original_dx', '')).strip().lower()
    if not original_dx or original_dx == 'nan':
        return 0

    ehr_primary_indicators = [
        'primary_malignant_neoplasm_of',
        'primary_malignant_neuroendocrine_neoplasm_of',
        'primary_basal_cell_carcinoma_of',
        'primary_low_grade_serous_adenocarcinoma_of',
        'primary_adenocarcinoma_of',
        'primary_squamous_cell_carcinoma_of',
        'primary_central_nervous_system',
        'primary_cutaneous',
        'primary_mediastinal_thymic_large_bcell_lymphoma',
        'fibrohistiocytic_tumor',
        'gastrointestinal_stromal_tumor',
        'infiltrating_ductal_carcinoma_of',
        'infiltrating_duct_carcinoma_of',
        'inflammatory_carcinoma_of',
        'cholangiocarcinoma_of_biliary_tract',
        'intrahepatic_bile_duct_carcinoma',
        'clear_cell_carcinoma_of_kidney',
        'renal_cell_carcinoma',
        'clear_cell_adenocarcinoma_of_ovary',
        'serous_cystadenoma_borderline_malignancy_of_ovary',
        'endometrial_intraepithelial_neoplasia',
        'acute_leukemia',
        'acute_lymphoid_leukemia',
        'acute_myeloid_leukemia',
        'acute_promyelocytic_leukemia',
        'chronic_leukemia',
        'chronic_lymphoid_leukemia',
        'chronic_myeloid_leukemia',
        'chronic_myelomonocytic_leukemia',
        'multiple_myeloma',
        'hodgkin_lymphoma',
        'follicular_low_grade_bcell_lymphoma',
        'mantle_cell_lymphoma',
        'marginal_zone_lymphoma',
        'malignant_melanoma_of',
        'microcystic_adnexal_carcinoma',
        'sebaceous_adenocarcinoma_of',
        'carcinoid_tumor_of',
        'medullary_thyroid_carcinoma',
        'follicular_thyroid_carcinoma',
        'papillary_thyroid_carcinoma',
        'squamous_cell_carcinoma_of',
        'carcinomainsitu_of',
        'cervical_intraepithelial_neoplasia',
        'liver_cell_carcinoma',
        'vaginal_intraepithelial_neoplasia',
        'vulval_intraepithelial_neoplasia',
        'burkitts_lymphoma',
        'adenocarcinoma_of_',
        'mesothelioma_malignant',
        'malignant_mesothelioma_of'
    ]
    
    # Split by semicolon and keep only terms containing "primary"
    diagnoses = [
        dx.strip() 
        for dx in original_dx.split(';') 
        if any(ind in dx.lower() for ind in ehr_primary_indicators)
    ]


    # Map each dx to YAML system(s), then flatten
    mapped_systems = []
    for dx in diagnoses:
        if dx in dx_to_systems:
            mapped_systems.extend(dx_to_systems[dx])

    # Count how many are in the PGC list
    return sum(sys.lower() in [c.lower() for c in pgc_cancers] for sys in mapped_systems)


# Apply to dataframe
df['Possibly_Genetic_Cancers'] = df.apply(count_pgc, axis=1)

# # --- Function to count PGC cancers per patient ---
# def count_pgc(row):
#     if row.get('cancer', '').strip().lower() == 'control':
#         return 0

#     grouped = str(row.get('original_dx_grouped', '')).strip().lower()
#     if not grouped or grouped == 'nan':
#         return 0

#     # Split by semicolon and count overlaps (case-insensitive)
#     diagnoses = [dx.strip().lower() for dx in grouped.split(';')]
#     return sum(dx in [c.lower() for c in pgc_cancers] for dx in diagnoses)

# # --- Apply function ---
# df['Possibly_Genetic_Cancers'] = df.apply(count_pgc, axis=1)

# --- Save result ---
df.to_csv('dfci-ufc.aou.phenos.v2.tsv', sep='\t', index=False)
