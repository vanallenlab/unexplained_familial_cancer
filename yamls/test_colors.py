import yaml
import matplotlib.pyplot as plt

# Load YAML file
with open("color_scheme.yaml") as f:
    color_dict = yaml.safe_load(f)

# Remove double-hash if present
for k, v in color_dict.items():
    color_dict[k] = v.replace("##", "#")

# Prepare figure
fig, ax = plt.subplots(figsize=(6, len(color_dict) * 0.4))
y_pos = list(range(len(color_dict)))
colors = list(color_dict.values())
labels = list(color_dict.keys())

# Plot horizontal bars
ax.barh(y_pos, [1]*len(color_dict), color=colors)
ax.set_yticks(y_pos)
ax.set_yticklabels(labels)
ax.set_xticks([])
ax.invert_yaxis()  # top to bottom

plt.title("Cancer Type Color Scheme")
plt.tight_layout()
plt.show()
