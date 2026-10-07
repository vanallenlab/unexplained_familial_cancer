#!/usr/bin/env python3
import argparse
import os

import hail as hl


def initialize_hail(google_project_id: str, ref_genome: str) -> None:
    hl.init(
        gcs_requester_pays_configuration=google_project_id
    )
    hl.default_reference(new_default_reference = ref_genome)


def filter_matrix_table(hail_mt: hl.MatrixTable, interval: str) -> hl.MatrixTable:
    mt_filtered = hl.filter_intervals(
        hail_mt,
        [hl.parse_locus_interval(interval)]
    )

    return mt_filtered


def main():
    parser = argparse.ArgumentParser(description="Filter Hail MatrixTable on GCP")
    parser.add_argument("--input", default="gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/exome/splitMT/hail.mt" required=True, help="GCS input MT path")
    parser.add_argument("--project", default=os.getenv("GOOGLE_CLOUD_PROJECT"), help="GCP project for requester-pays access")
    parser.add_argument("--genome", default="GRCh38", help="Reference genome")
    parser.add_argument("--intervals-file")
    args = parser.parse_args()

   

    print(f"[filter_intervals] starting: input={args.input}", flush=True)

    initialize_hail(google_project_id=args.project, ref_genome="GRCh38")
    print("[filter_intervals] hail initialized", flush=True)

    mt = hl.read_matrix_table(args.input)

    with open(args.intervals_file) as f:
        for line in f:
            chromosome, start, end, gene = line.rstrip("\n").split("\t")
            interval = chromosome + ":" + start + "-" + end
            filtered_mt = filter_matrix_table(hail_mt=mt, interval=interval)
            hl.export_vcf(filtered_mt, f"gs://dataproc-staging-wb-cordial-diamond-9893/genes/{gene}.vcf.bgz")

    print("[filter_intervals] interval filter applied", flush=True)


if __name__ == "__main__":
    main()
