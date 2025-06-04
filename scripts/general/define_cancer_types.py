import pandas as pd
import yaml
from collections import defaultdict

# Load data
tsv_path = "dfci-ufc.aou.phenos.tsv"
yaml_path = "cancer_ontology.yaml"
df = pd.read_csv(tsv_path, sep="\t")
df["original_dx"] = df["original_dx"].fillna("")

# Parse YAML ontology
with open(yaml_path, "r") as f:
    ontology = yaml.safe_load(f)

# Maps each diagnosis string to the groups it belongs to
string_to_groups = defaultdict(set)

# Maps each group to all diagnosis strings under it (direct + nested)
group_to_strings = defaultdict(set)

def collect_types(system, parent_groups=[]):
    current_group = system.get("name")
    full_group_hierarchy = parent_groups + [current_group]

    # Add each dx string to all parent groups
    for dx in system.get("types", []):
        dx = dx.strip()
        for group in full_group_hierarchy:
            string_to_groups[dx].add(group)
            group_to_strings[group].add(dx)

    # Recurse into sub-systems
    for sub in system.get("sub_systems", []):
        collect_types(sub, full_group_hierarchy)

# Traverse full ontology
for system in ontology["cancer_ontology"]["systems"]:
    collect_types(system)

# Track patients with each dx and each group
string_counts = defaultdict(set)
group_counts = defaultdict(set)

# Count patients per diagnosis string and group
for _, row in df.iterrows():
    sample_id = row["Sample"]
    dx_list = [s.strip() for s in row["original_dx"].split(";") if s.strip()]
    for dx in dx_list:
        if dx in string_to_groups:
            string_counts[dx].add(sample_id)
            for group in string_to_groups[dx]:
                group_counts[group].add(sample_id)

def redact(count):
    return "<20" if count < 20 else str(count)

# Print results
print("\n🔹 Diagnosis-level patient counts:")
for dx, samples in sorted(string_counts.items()):
    print(f"{dx}: {redact(len(samples))}")

print("\n🔹 Group-level patient counts (including nested):")
for group, samples in sorted(group_counts.items()):
    print(f"{group}: {redact(len(samples))}")

