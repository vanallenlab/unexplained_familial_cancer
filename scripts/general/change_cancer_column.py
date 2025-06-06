import yaml
import pandas as pd

# --- Load YAML file ---
with open('cancer_ontology.yaml', 'r') as f:
    systems = yaml.safe_load(f)

# --- Build a flat mapping from diagnosis to system labels ---
dx_to_systems = {}

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
            systems_found.add(';'.join(dx_to_systems[dx]))

    if systems_found:
        return '|'.join(sorted(systems_found))
    else:
        return 'NA'

# --- Apply function to create new column ---
df['original_dx_grouped'] = df.apply(group_dxs, axis=1)

# --- Save result ---
df.to_csv('dfci-ufc.aou.phenos.adjusted.tsv', sep='\t', index=False)
