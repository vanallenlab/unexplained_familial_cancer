# The Unexplained Familial Cancer (UFC)
# Copyright (c) 2023-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0g

version 1.0
#import "Ufc_utilities/Ufc_utilities.wdl" as Tasks

workflow REPLICATION_1_EXTRACT_GENOMIC_RANGES {
	input {
		String AllofUs_HailMatrix = "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/exome/splitMT/hail.mt"
		String Google_Project = "wb-cordial-diamond-9893"
		String AllofUs_version = "CDRv9"
		File MANE_GENCODE

		File Gene_List
		Int Gene_Buffer

		String project_name
	}

	Array[String] Genes = read_lines(Gene_List)
	
	call T1_Find_Gene_Ranges {
		input:
			MANE_GENCODE = MANE_GENCODE,
			project_name = project_name,
			Genes = Genes,
			Gene_Buffer = Gene_Buffer
	}

	scatter (gene_range in T1_Find_Gene_Ranges.out1) {

		String interval = gene_range[0]
		String gene = gene_range[1]

		call T2_Query_Hail_MT {
			input:
				AllofUs_HailMatrix = AllofUs_HailMatrix,
				Gene = gene,
				Interval = interval,
				project_name = project_name,
				Google_Project = Google_Project
		}

		call T3_Write_Index {
			input:
				vcf = T2_Query_Hail_MT.out1
		}
	}
}

# This task takes genes and finds their 
# associated genomic intervals which can 
# be used to later extract in Hail MT

task T1_Find_Gene_Ranges {
	input {
		File MANE_GENCODE
		String project_name

		Array[String] Genes
		Int Gene_Buffer
	}
	String output_file = "~{project_name}.genomic_intervals.bed"

	command <<<
	set -euxo pipefail

	echo "STARTING GUNZIP" >&2
	gunzip -c ~{MANE_GENCODE} > MANE.GRCh38.gtf
	echo "GUNZIP FINISHED" >&2

	echo "STARTING PYTHON" >&2
	python3 <<CODE
	import pandas as pd

	f = open("~{output_file}","w")

	mane_df = pd.read_csv("MANE.GRCh38.gtf",
		sep='\t',
		index_col=False,
		comment='#',
		header=None,
		names = ["chrom","source","feature","start","end","score","strand","frame","attributes"])

	genes = "~{sep=' ' Genes}".split()

	print("hi")
	# Loop through genes
	for gene in genes:
		gene_df = mane_df[mane_df["attributes"].str.contains(f"gene_name \"{gene}\"", na=False)]
		gene_df = gene_df[gene_df['feature'] == "gene"]

		if gene_df.empty:
			raise ValueError(f"Could not find gene: {gene}")

		row = gene_df.iloc[0]
		f.write((row['chrom']) + ":" + str(int(row['start']) - ~{Gene_Buffer}) + "-" + str(int(row['end']) + ~{Gene_Buffer}) + "\t" + gene + "\n")


	CODE
	>>>
	runtime {
		docker: "vanallenlab/hail_rf"
		memory: "8 GB"
	}
	output {
		Array[Array[String]] out1 = read_tsv(output_file)
	}
}

task T2_Query_Hail_MT {
	input {
		String AllofUs_HailMatrix
		String project_name
		String Google_Project

		String Gene
		String Interval
	}

	command <<<
	set -euxo pipefail

	python3 <<CODE
	import hail as hl
	import os

	#print(os.getenv())
	for key, value in sorted(os.environ.items()):
		print(f"{key}={value}")

	hl.init(gcs_requester_pays_configuration="~{Google_Project}")
	hl.default_reference(new_default_reference = "GRCh38")

	mt = hl.read_matrix_table("~{AllofUs_HailMatrix}")
	mt_filtered = hl.filter_intervals(mt,[hl.parse_locus_interval("~{Interval}")])

	# Keep only variant rows
	#mt_filtered = mt_filtered.select_entries()
	variant_table = mt_filtered.rows()

	#hl.export_vcf(mt_filtered,f"gs://dataproc-staging-wb-cordial-diamond-9893/genes/~{Gene}.vcf.bgz")
	hl.export_vcf(variant_table,"~{Gene}.vcf.bgz")
	CODE
	>>>
	runtime {
		docker: "vanallenlab/hail_rf"
		memory: "8 GB"
		preemptible: 3
	}
	output {
		File out1 = "~{Gene}.vcf.bgz"
	}
}

task T3_Write_Index {
	input {
		File vcf
	}

	String rootname = basename(vcf)

	command <<<
	set -euxo pipefail

	echo "Indexing ~{rootname}"

	bcftools index -t "~{vcf}"

	echo "Uploading index"

	export CLOUDSDK_PYTHON=/opt/conda/envs/gatk-sv/bin/python3

	gcloud storage cp ~{vcf} gs://dataproc-staging-wb-cordial-diamond-9893/genes_raw_vcf/
	gcloud storage cp ~{vcf}.tbi gs://dataproc-staging-wb-cordial-diamond-9893/genes_raw_vcf/

	echo "Finished ~{rootname}"
	>>>

	runtime {
		docker: "vanallenlab/g2c_pipeline"
		memory: "8 GB"
		preemptible: 3
	}
}
