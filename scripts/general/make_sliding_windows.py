import sys

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

BUFFER_BP = 2_000_000  # 2 Mb buffer on each side


# =========================================================
# Load & merge masked regions (new format)
# =========================================================
def load_buffered_regions(region_file):
    """
    Reads tab-delimited regions:
    col2 = chrom, col3 = start, col4 = end
    Applies ±2 Mb buffer and merges overlaps.
    Returns: dict chrom -> list of (start, end)
    """
    regions = {}

    with open(region_file) as f:
        for line in f:
            if not line.strip():
                continue

            fields = line.rstrip().split("\t")
            chrom = fields[1]
            start = int(fields[2])
            end = int(fields[3])

            if chrom not in CHROM_SIZES:
                continue

            chrom_len = CHROM_SIZES[chrom]

            buffered_start = max(1, start - BUFFER_BP)
            buffered_end = min(chrom_len, end + BUFFER_BP)

            regions.setdefault(chrom, []).append((buffered_start, buffered_end))

    # Merge overlapping buffered intervals per chromosome
    merged = {}

    for chrom, intervals in regions.items():
        intervals.sort()
        merged_intervals = []

        cur_start, cur_end = intervals[0]

        for s, e in intervals[1:]:
            if s <= cur_end:
                cur_end = max(cur_end, e)
            else:
                merged_intervals.append((cur_start, cur_end))
                cur_start, cur_end = s, e

        merged_intervals.append((cur_start, cur_end))
        merged[chrom] = merged_intervals

    return merged


# =========================================================
# Check overlap between window and masked regions
# =========================================================
def overlaps_any(win_start, win_end, intervals):
    for s, e in intervals:
        if win_start <= e and win_end >= s:
            return True
    return False


# =========================================================
# Main sliding window generator
# =========================================================
def generate_sliding_windows(
    window_size_bp,
    step_size_bp,
    region_file,
    output_filename="genome_windows.txt"
):
    masked = load_buffered_regions(region_file)

    total_kept = 0
    total_independent = 0

    with open(output_filename, "w") as out:
        print(f"Generating windows (window={window_size_bp}, step={step_size_bp})")
        print("Excluding windows overlapping buffered regions (±2 Mb)\n")

        sorted_chrs = sorted(CHROM_SIZES.keys(), key=lambda x: (len(x), x))

        for chrom in sorted_chrs:
            chrom_len = CHROM_SIZES[chrom]
            mask_list = masked.get(chrom, [])

            start = 1
            kept = 0
            independent = 0

            while start + window_size_bp - 1 <= chrom_len:
                end = start + window_size_bp - 1
                window_id = f"{chrom}_{start}_{end}"

                if not overlaps_any(start, end, mask_list):
                    out.write(
                        f"{chrom}\t{start}\t{end}\t{window_id}\t{window_size_bp}\n"
                    )
                    kept += 1

                    # independent windows (non-overlapping grid)
                    if (start - 1) % window_size_bp == 0:
                        independent += 1

                start += step_size_bp

            total_kept += kept
            total_independent += independent
            print(f"{chrom}: {kept} windows kept ({independent} independent)")

    print("\n====================================================")
    print(f"TOTAL WINDOWS KEPT: {total_kept}")
    print(f"TOTAL INDEPENDENT (NON-OVERLAPPING): {total_independent}")
    print(f"BONFERRONI THRESHOLD = 0.05 / {total_independent:.0f}")
    print("====================================================")
    print(f"Output written to: {output_filename}")


# =========================================================
# Run if script called directly
# =========================================================
if __name__ == "__main__":
    WINDOW_SIZE = 100_000   # 100 kb
    STEP_SIZE = 50_000      # 50 kb
    REGION_FILE = "reference_files/centromeres.txt"

    generate_sliding_windows(WINDOW_SIZE, STEP_SIZE, REGION_FILE, output_filename=f"genome_windows.{WINDOW_SIZE}bp.txt")
