import pandas as pd
import numpy as np
import glob
import os

def round_sig(x, sig=3):
    if pd.isna(x) or x == 0:
        return x
    return float(f"{x:.{sig}g}")

#out_file = "roh_results_summary.xlsx"
out_file = "Supplementary Table 7.sw_100kb.roh_500kb.xlsx"

with pd.ExcelWriter(out_file, engine="openpyxl") as writer:

    for fn in sorted(glob.glob("*.tsv.gz")):
        print(fn)
        df = pd.read_csv(fn, sep="\t", compression="gzip")

        # keep only the first four columns
        df = df.iloc[:, :4]

        # round all columns except haplotype
        for col in df.columns:
            if col == "haplotype":
                continue

            if pd.api.types.is_numeric_dtype(df[col]):
                df[col] = df[col].apply(round_sig)

        sheet_name = os.path.basename(fn).replace(
            ".sw_250000.roh_250000.tsv.gz", ""
        )[:31]

        df.to_excel(writer, sheet_name=sheet_name, index=False)

print(f"Wrote {out_file}")
