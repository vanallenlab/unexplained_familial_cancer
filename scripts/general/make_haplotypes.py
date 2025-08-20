import pandas as pd

# Read the recombination map
df = pd.read_csv("genetic_map_hg38_withX.txt", delim_whitespace=True)

# Ensure correct column names
df.columns = df.columns.str.strip()
df = df.rename(columns={df.columns[0]: "chr", df.columns[1]: "position"})

# Function to compute rolling mean within ±range_bp
def local_mean(sub_df, col="COMBINED_rate(cM/Mb)", range_bp=1000):
    vals = []
    pos = sub_df["position"].values
    arr = sub_df[col].values

    for i, p in enumerate(pos):
        mask = (pos >= p - range_bp) & (pos <= p + range_bp)
        vals.append(arr[mask].mean())
    return pd.Series(vals, index=sub_df.index)

# Apply per chromosome
df["rolling_mean"] = (
    df.groupby("chr", group_keys=False)
      .apply(lambda g: local_mean(g, range_bp=1000))
)

# Example: define a breakpoint where rolling_mean > threshold
BREAK_THRESHOLD = 10  # or whatever threshold makes sense
df['is_break'] = df['rolling_mean'] > BREAK_THRESHOLD
df.to_csv("rolling_mean_genetic_map.txt",sep='\t',index=False)

# Read in precomputed rolling mean data
df = pd.read_csv("rolling_mean_genetic_map.txt", sep='\t', index_col=False)

def collapse_breaks_with_edges(subdf, min_len=10_000):
    """
    Collapse consecutive is_break=True rows into representative haplotypes,
    handle first/last haplotype based on non-zero COMBINED_rate(cM/Mb),
    and only keep haplotypes >= min_len.
    """
    if subdf.empty:
        return pd.DataFrame()

    filtered_rows = []
    block_rows = []

    # Identify positions with non-zero rate
    nonzero_positions = subdf[subdf['COMBINED_rate(cM/Mb)'] != 0]['position'].values
    if len(nonzero_positions) == 0:
        return pd.DataFrame()  # skip if entire chromosome is zero

    chrom_start = nonzero_positions[0]
    chrom_end = nonzero_positions[-1]

    for _, row in subdf.iterrows():
        pos = row['position']
        if row['is_break'] and chrom_start <= pos <= chrom_end:
            block_rows.append(row)
        else:
            if block_rows:
                best_row = max(block_rows, key=lambda r: r['COMBINED_rate(cM/Mb)'])
                filtered_rows.append(best_row)
                block_rows = []

    if block_rows:
        best_row = max(block_rows, key=lambda r: r['COMBINED_rate(cM/Mb)'])
        filtered_rows.append(best_row)

    # Create BED-like intervals
    bed_rows = []
    prev_end = chrom_start
    for r in filtered_rows:
        start = prev_end
        end = r['position']
        if end - start >= min_len:  # Only keep if >= min_len
            bed_rows.append({
                'chr': r['chr'],
                'start': start,
                'end': end,
                'name': f"{r['chr']}_{start}_{end}",
                'length': end - start
            })
        prev_end = end

    # Handle final haplotype to chromosome end
    if prev_end < chrom_end and (chrom_end - prev_end >= min_len):
        bed_rows.append({
            'chr': subdf.iloc[0]['chr'],
            'start': prev_end,
            'end': chrom_end,
            'name': f"{subdf.iloc[0]['chr']}_{prev_end}_{chrom_end}",
            'length': chrom_end - prev_end
        })

    return pd.DataFrame(bed_rows)

# Apply per chromosome
collapsed_bed = []
for chrom, subdf in df.groupby('chr'):
    collapsed_bed.append(collapse_breaks_with_edges(subdf, min_len=10_000))

bed_df = pd.concat(collapsed_bed).sort_values(['chr', 'start'])

# Write BED file with length column
bed_df.to_csv("collapsed_haplotypes_min10kb.bed", sep='\t', index=False, header=False)

