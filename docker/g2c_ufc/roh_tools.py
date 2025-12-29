import pandas as pd
import numpy as np
import argparse

# =========================================================
# Chromosome sizes (GRCh38)
# =========================================================
CHROM_SIZES = {
    'chr1': 248956422,
    'chr2': 242193529,
    'chr3': 198295559,
    'chr4': 190214555,
    'chr5': 181538259,
    'chr6': 170805979,
    'chr7': 159345973,
    'chr8': 145138636,
    'chr9': 138394717,
    'chr10': 133797422,
    'chr11': 135086622,
    'chr12': 133275309,
    'chr13': 114364328,
    'chr14': 107043718,
    'chr15': 101991189,
    'chr16': 90338345,
    'chr17': 83257441,
    'chr18': 80373285,
    'chr19': 58617616,
    'chr20': 64444167,
    'chr21': 46709983,
    'chr22': 50818468,
}
TOTAL_AUTOSOMAL_SIZE = sum(CHROM_SIZES.values())
roh_thresholds = ['1000','10000','50000','100000','250000','500000','1000000','3000000','5000000']
def build_roh_fraction_table(roh_df,roh_thresholds):
	# -----------------------
	# Aggregate ROH length per sample per threshold
	# -----------------------
	out = pd.DataFrame({"Sample": roh_df["Sample"].unique()}).set_index("Sample")

	for thresh in roh_thresholds:
	    colname = f"ROH_sum_ge_{thresh}"

	    summed = (
	        roh_df[roh_df["Length"] >= thresh]
	        .groupby("Sample")["Length"]
	        .sum()
	    )

	    out[colname] = summed

	for thresh in roh_thresholds:
	    col_name = f"ROH_frac_ge_{thresh}"

	    roh_sum = (
	        roh_df[roh_df["Length"] >= thresh]
	        .groupby("Sample")["Length"]
	        .sum()
	    )

	    out[col_name] = roh_sum / TOTAL_AUTOSOMAL_SIZE

	    

	# Replace NaNs with 0 (no ROH meeting threshold)
	out = out.fillna(0).reset_index()

	# -----------------------
	# Save
	# -----------------------
	out.to_csv("roh_stats/percent_roh_per_person.tsv", sep="\t", index=False)
	return out

# -----------------------
# Function 2:
# Phenotype-aware summaries + stats
# -----------------------
def roh_phenotype_analysis(roh_frac_df, phenotype_df):
    phenotype_df["Sample"] = phenotype_df["Sample"].astype(str)

    df = roh_frac_df.merge(phenotype_df, on="Sample", how="left")

    results = []

    # ---- Excess ROH (>3% at ≥1MB) ----
    excess_mask = df["ROH_frac_ge_1000000"] > 0.03
    excess_ids = df.loc[excess_mask, "Sample"].tolist()
    excess_pct = excess_mask.mean() * 100

    results.append({
        "metric": "excess_ROH_>3pct_ge_1MB",
        "value": f"{excess_pct:.2f}%",
        "details": ",".join(excess_ids)
    })

    # ---- Mean ROH fraction at each threshold ----
    for thresh in roh_thresholds:
        col = f"ROH_frac_ge_{thresh}"
        results.append({
            "metric": f"mean_{col}",
            "value": df[col].mean(),
            "details": ""
        })

    summary_df = pd.DataFrame(results)
    summary_df.to_csv("roh_summary_metrics.tsv", sep="\t", index=False)

    # ---- Log file by ancestry ----
    log_rows = []
    for thresh in ROH_THRESHOLDS:
        col = f"ROH_frac_ge_{thresh}"
        for pop, g in df.groupby("intake_qc_pop"):
            log_rows.append({
                "threshold": thresh,
                "ancestry": pop,
                "mean_frac": g[col].mean(),
                "n": len(g)
            })

    log_df = pd.DataFrame(log_rows)
    log_df.to_csv("roh_stats/roh_by_ancestry_log.tsv", sep="\t", index=False)

    # ---- ANOVA by ancestry ----
    anova_rows = []
    for thresh in ROH_THRESHOLDS:
        col = f"ROH_frac_ge_{thresh}"

        groups = [
            g[col].values
            for _, g in df.groupby("intake_qc_pop")
            if len(g) >= 5
        ]

        if len(groups) >= 2:
            fstat, pval = stats.f_oneway(*groups)
            anova_rows.append({
                "threshold": thresh,
                "anova_p": pval,
                "anova_F": fstat
            })

    anova_df = pd.DataFrame(anova_rows)
    anova_df.to_csv("roh_stats/roh_anova_by_ancestry.tsv", sep="\t", index=False)

# -----------------------
# Function 3:
# ROH burden characteristics
# -----------------------
def roh_burden_stats(roh_df):
    rows = []

    for thresh in roh_thresholds:
        sub = roh_df[roh_df["Length"] >= thresh]

        per_sample = sub.groupby("Sample").size()
        lengths = sub["Length"]

        rows.append({
            "threshold": thresh,
            "mean_roh_per_sample": per_sample.mean(),
            "sd_roh_per_sample": per_sample.std(),
            "mean_roh_length": lengths.mean(),
            "sd_roh_length": lengths.std(),
            "total_rohs": len(sub)
        })

    out = pd.DataFrame(rows)
    out.to_csv("roh_stats/roh_burden_stats.tsv", sep="\t", index=False)

def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--phenotype-data', required=True, help='data describing phenotypes')
    parser.add_argument('--samples', required=True,help='list of samples to include')
    parser.add_argument('--roh-file', required=True, help='bcftools roh filtered output')

    args = parser.parse_args()

    # Load ROH data and phenotypic data
    roh_df = pd.read_csv(args.roh_file,sep='\t',index=False,header=None,names=['RG','Sample','Chromosome','Start','End','Length','Num_markers','Quality'])
    phenotype_df = pd.read_csv("phenotype_data.tsv",sep="\t", dtype=str)

    # Ensure correct dtypes
    roh_df["Sample"] = roh_df["Sample"].astype(str)
	roh_df["Length"] = roh_df["Length"].astype(float)

    # Load samples list (one sample ID per line)
    samples = pd.read_csv(args.samples,header=None,dtype=str)[0].tolist()

    samples = set(samples)

    # Filter ROH dataframe
    roh_df = roh_df[roh_df["Sample"].isin(samples)]

    

    # Run everything
    roh_frac_df = build_roh_fraction_table(roh_df)
    roh_phenotype_analysis(roh_frac_df, phenotype_df)
    roh_burden_stats(roh_df)




if __name__ == '__main__':
  main()