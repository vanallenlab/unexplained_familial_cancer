import pandas as pd
import matplotlib.pyplot as plt

# ------------------------------------------------------------
# Load ROH per-base signal
# ------------------------------------------------------------
df = pd.read_csv("roh_per_base_normalized_signal.tsv.gz", sep="\t")

# ------------------------------------------------------------
# Smooth signals (rolling mean)
# Window size depends on resolution; adjust as needed
# ------------------------------------------------------------
window = 500  # e.g. 500 bp smoothing

df["case_smooth"] = df["case_roh_fraction"].rolling(
    window=window,
    center=True,
    min_periods=1
).mean()

df["control_smooth"] = df["control_roh_fraction"].rolling(
    window=window,
    center=True,
    min_periods=1
).mean()

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------
plt.figure(figsize=(12, 4))

plt.plot(
    df["pos"],
    df["case_smooth"],
    label="Cases",
    linewidth=2,
    alpha=0.9
)

plt.plot(
    df["pos"],
    df["control_smooth"],
    label="Controls",
    linewidth=2,
    alpha=0.9
)

plt.xlabel("Genomic position (bp)")
plt.ylabel("Fraction of samples with overlapping ROH")
plt.title("Normalized ROH burden across haplotype")

plt.legend()
plt.tight_layout()
plt.show()
