
# The Unexplained Familial Cancer (UFC)
# Copyright (c) 2023-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

version 1.0

workflow EXTRACT_ACAF_RANGES_AoU {
    input {
        # Must point to an ACAF-filtered Hail MatrixTable.
        String AllofUs_HailMatrix

        String Google_Project = "wb-cordial-diamond-9893"
        String AllofUs_version = "CDRv9"

        String region_name

        # Supported values: is_snv, is_region
        String Analysis_Type

        # SNV:    chr1_3763798_C_T
        # Region: chr1:3763798-3768798
        String Target

        Int Gene_Buffer = 1000000

        Float LD_R2_Threshold = 0.2
    }

    call T1_Query_Hail_MT {
        input:
            AllofUs_HailMatrix = AllofUs_HailMatrix,
            Google_Project = Google_Project,
            Target = Target,
            Analysis_Type = Analysis_Type,
            Gene_Buffer = Gene_Buffer,
            region_name = region_name
    }

    call T2_Write_Index {
        input:
            vcf = T1_Query_Hail_MT.out1,
            region_name = region_name
    }

    if (Analysis_Type == "is_snv") {
        call T3_PLINK_LD {
            input:
                vcf = T1_Query_Hail_MT.out1,
                Target = Target,
                region_name = region_name,
                LD_R2_Threshold = LD_R2_Threshold
        }
    }

}


# Parse the target, expand the interval, validate SNVs,
# query the ACAF MatrixTable, and export the genotype data.
task T1_Query_Hail_MT {
    input {
        String AllofUs_HailMatrix
        String Google_Project
        String Target
        String Analysis_Type
        Int Gene_Buffer
        String region_name
    }

    command <<<
    set -euo pipefail

    python3 << CODE
    import re
    import hail as hl

    target = "~{Target}"
    mode = "~{Analysis_Type}"
    buffer = ~{Gene_Buffer}

    if buffer < 0:
        raise ValueError("Gene_Buffer must be nonnegative")

    # Parse target and expand by the requested buffer.
    if mode == "is_snv":
        match = re.fullmatch(
            r"chr([0-9]+|X|Y)_([0-9]+)_([ACGT]+)_([ACGT]+)",
            target
        )
        if not match:
            raise ValueError(
                f"Invalid SNV identifier: {target}. "
                "Expected chr1_3763798_C_T."
            )

        chrom, pos, ref, alt = match.groups()
        pos = int(pos)

        if ref == alt:
            raise ValueError("Reference and alternate alleles must differ")

        start = max(1, pos - buffer)
        end = pos + buffer

    elif mode == "is_region":
        match = re.fullmatch(
            r"chr([0-9]+|X|Y):([0-9]+)-([0-9]+)",
            target
        )
        if not match:
            raise ValueError(
                f"Invalid region: {target}. "
                "Expected chr1:3763798-3768798."
            )

        chrom, start, end = match.groups()
        start, end = int(start), int(end)

        if start < 1 or end < start:
            raise ValueError(f"Invalid genomic coordinates: {target}")

        start = max(1, start - buffer)
        end = end + buffer

    else:
        raise ValueError(
            f"Unsupported Analysis_Type: {mode}. "
            "Use is_snv or is_region."
        )

    interval = f"{chrom}:{start}-{end}"

    print(f"Target: {target}")
    print(f"Analysis type: {mode}")
    print(f"Extracting GRCh38 interval: {interval}")

    hl.init(
        gcs_requester_pays_configuration="~{Google_Project}"
    )
    hl.default_reference("GRCh38")

    mt = hl.read_matrix_table("~{AllofUs_HailMatrix}")

    hl_interval = hl.parse_locus_interval(
        interval,
        reference_genome="GRCh38"
    )

    mt_filtered = hl.filter_intervals(mt, [hl_interval])

    # Confirm that the requested SNV exists with the expected alleles.
    if mode == "is_snv":
        rows = mt_filtered.rows()

        target_locus = hl.parse_locus(
            f"{chrom}:{pos}",
            reference_genome="GRCh38"
        )

        target_rows = rows.filter(
            (rows.locus == target_locus) &
            (rows.alleles[0] == ref) &
            (rows.alleles[1] == alt)
        )

        n_target = target_rows.count()

        if n_target == 0:
            raise ValueError(
                f"Target SNV not found in the ACAF MatrixTable: {target}"
            )

        print(f"Validated target SNV: {target}")

    # Export the MatrixTable to preserve genotype entries.
    output_vcf = "~{region_name}.vcf.bgz"
    hl.export_vcf(mt_filtered, output_vcf)

    with open("genomic_interval.txt", "w") as f:
        f.write(interval + "\n")

    print(f"Exported: {output_vcf}")
    CODE
    >>>

    output {
        File out1 = "~{region_name}.vcf.bgz"
        String Interval = read_string("genomic_interval.txt")
    }

    runtime {
        docker: "vanallenlab/hail_rf"
        memory: "8 GB"
        preemptible: 3
    }
}


# Index the extracted VCF and upload the VCF and index.
task T2_Write_Index {
    input {
        File vcf
        String region_name
    }

    command <<<
    set -euo pipefail

    echo "Indexing ~{vcf}"
    bcftools index -t "~{vcf}"
    >>>

    output {
        File index = "~{vcf}.tbi"
    }

    runtime {
        docker: "vanallenlab/g2c_pipeline"
        preemptible: 3
    }
}


# For SNV targets, identify variants in LD with the target.
task T3_PLINK_LD {
    input {
        File vcf
        String Target
        String region_name
        Float LD_R2_Threshold
    }

    command <<<
    set -euxo pipefail

    python3 << CODE
    import re

    target = "~{Target}"

    chrom, pos, ref, alt = target.split('_')

    # PLINK variant ID convention: CHROM:POS:REF:ALT
    target_id = f"{chrom}:{pos}:{ref}:{alt}"

    with open("target_id.txt", "w") as f:
        f.write(target_id + "\n")

    CODE

    target_id=$(cat target_id.txt)

    plink \
        --vcf "~{vcf}" \
        --double-id \
        --allow-extra-chr \
        --set-missing-var-ids '@:#:$r:$a' \
        --ld-snp "$target_id" \
        --r2 \
        --ld-window-kb 1000 \
        --ld-window 99999 \
        --ld-window-r2 ~{LD_R2_Threshold} \
        --out "~{region_name}.ld"

    test -s "~{region_name}.ld.ld"

    echo "LD analysis completed for ${target_id}"
    >>>

    output {
        File ld_results = "~{region_name}.ld.ld"
    }

    runtime {
        docker: "vanallenlab/g2c_pipeline"
        memory: "8 GB"
        preemptible: 3
    }
}
