import plotly.graph_objects as go

labels = [
    "Malignant Neoplastic Disease (49,181)",
    "srWGS (34,932)",
    "Excluded before srWGS",
    "Responded to FDR cancer data (11,933)",
    "No FDR response",
    "Female (breast only)",
    "Breast diagnosis",
    "Kidney diagnosis",
    "Lung diagnosis",
    "Breast: 3 FDRs with cancer",
    "Kidney: 3 FDRs with cancer",
    "Lung: 3 FDRs with cancer"
]

idx = {name: i for i, name in enumerate(labels)}

source = []
target = []
value = []

# 49,181 → srWGS + excluded
source += [idx["Malignant Neoplastic Disease (49,181)"]]
target += [idx["srWGS (34,932)"]]
value  += [34932]

source += [idx["Malignant Neoplastic Disease (49,181)"]]
target += [idx["Excluded before srWGS"]]
value  += [49181 - 34932]

# 34,932 → FDR responders + non-responders
source += [idx["srWGS (34,932)"]]
target += [idx["Responded to FDR cancer data (11,933)"]]
value  += [11933]

source += [idx["srWGS (34,932)"]]
target += [idx["No FDR response"]]
value  += [34932 - 11933]

# breast path
source += [idx["Responded to FDR cancer data (11,933)"]]
target += [idx["Female (breast only)"]]
value  += [7083]

source += [idx["Female (breast only)"]]
target += [idx["Breast diagnosis"]]
value  += [2439]

source += [idx["Breast diagnosis"]]
target += [idx["Breast: 3 FDRs with cancer"]]
value  += [390]

# kidney path
source += [idx["Responded to FDR cancer data (11,933)"]]
target += [idx["Kidney diagnosis"]]
value  += [391]

source += [idx["Kidney diagnosis"]]
target += [idx["Kidney: 3 FDRs with cancer"]]
value  += [40]

# lung path
source += [idx["Responded to FDR cancer data (11,933)"]]
target += [idx["Lung diagnosis"]]
value  += [987]

source += [idx["Lung diagnosis"]]
target += [idx["Lung: 3 FDRs with cancer"]]
value  += [87]

fig = go.Figure(go.Sankey(
    node=dict(
        pad=20,
        thickness=20,
        label=labels
    ),
    link=dict(
        source=source,
        target=target,
        value=value
    )
))

fig.update_layout(
    title="Cancer Cohort Filtering Flow",
    font_size=12
)

fig.show()