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


# =========================================================
# Load & merge N-masked regions
# =========================================================
def load_nmasked(mask_file, min_required_length=50000):
    """
    Reads N-masked regions and merges overlapping or adjacent intervals.
    Keeps only intervals with length >= min_required_length.
    Returns: dict of chrom → list of (start, end) merged intervals
    """
    masked = {}

    with open(mask_file) as f:
        for line in f:
            fields = line.strip().split("\t")
            if len(fields) < 7:
                continue

            chrom = fields[1]
            start = int(fields[2])
            end = int(fields[3])
            length = int(fields[6])

            # keep only long N-mask intervals
            if length < min_required_length:
                continue

            masked.setdefault(chrom, []).append((start, end))

    # Merge each chromosome's intervals
    merged = {}

    for chrom, intervals in masked.items():
        if not intervals:
            continue

        intervals.sort()
        merged_intervals = []

        current_start, current_end = intervals[0]

        for s, e in intervals[1:]:
            if s <= current_end:  # overlapping or adjacent
                current_end = max(current_end, e)
            else:
                merged_intervals.append((current_start, current_end))
                current_start, current_end = s, e

        merged_intervals.append((current_start, current_end))
        merged[chrom] = merged_intervals

    return merged


# =========================================================
# Compute overlap between a window and N-masked intervals
# =========================================================
def masked_overlap(win_start, win_end, intervals):
    total = 0
    for s, e in intervals:
        overlap = max(0, min(win_end, e) - max(win_start, s) + 1)
        total += overlap
    return total


# =========================================================
# Main sliding window generator
# =========================================================
def generate_sliding_windows(window_size_bp, step_size_bp, mask_file, output_filename="genome_windows.txt"):
    nmasked = load_nmasked(mask_file, min_required_length=50000)

    total_independent = 0
    total_windows_kept = 0

    with open(output_filename, "w") as f:
        print(f"Generating Windows (window={window_size_bp}, step={step_size_bp})")
        print(f"Excluding windows with ≥50kb overlap with N-masked regions\n")

        sorted_chrs = sorted(CHROM_SIZES.keys(), key=lambda x: (len(x), x))

        for chrom in sorted_chrs:
            chrom_length = CHROM_SIZES[chrom]
            mask_list = nmasked.get(chrom, [])

            start = 1
            kept = 0
            independent = 0

            while start + window_size_bp - 1 <= chrom_length:
                end = start + window_size_bp - 1
                window_id = f"{chrom}_{start}_{end}"
                length = end - start + 1

                # compute overlap with N-masked regions
                overlap = masked_overlap(start, end, mask_list)

                if overlap < 50000:
                    f.write(f"{chrom}\t{start}\t{end}\t{window_id}\t{length}\n")
                    kept += 1

                    # Non-overlapping windows:
                    # windows starting at positions 1, 1+window, 1+2*window, ...
                    if (start - 1) % window_size_bp == 0:
                        independent += 1

                start += step_size_bp

            total_windows_kept += kept
            total_independent += independent

            print(f"{chrom}: {kept} windows kept ({independent} independent)")

    print("\n====================================================")
    print(f"TOTAL WINDOWS KEPT: {total_windows_kept}")
    print(f"TOTAL INDEPENDENT (NON-OVERLAPPING): {total_independent}")
    print(f"BONFERRONI THRESHOLD = 0.05 / {total_independent:.0f}")
    print("====================================================")
    print(f"Output written to: {output_filename}")


# =========================================================
# Run if script called directly
# =========================================================
if __name__ == "__main__":
    WINDOW_SIZE = 100000    # 100 kb
    STEP_SIZE = 50000       # 50 kb
    MASK_FILE = "gap.txt"

    generate_sliding_windows(WINDOW_SIZE, STEP_SIZE, MASK_FILE)
