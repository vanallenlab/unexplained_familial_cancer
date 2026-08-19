# The Unexplained Familial Cancer (UFC)
# Copyright (c) 2023-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0g

version 1.0
#import "Ufc_utilities/Ufc_utilities.wdl" as Tasks

workflow REPLICATION_1_EXTRACT_GENOMIC_RANGES {
	input {
		String AllofUs_HailMatrix
		String Google_Bucket_ID
		String AllofUs_version = "CDRv9"
		File MANE_GENCODE

		Array[String] Genes
		Int Gene_Buffer
		Array[String]? SNPs
		Array[String]? SNP_names
		String? SNP_buffer

		String project_name
	}

	call T1_Find_Gene_Ranges {
		input:
	  		MANE_GENCODE = MANE_GENCODE,
	  		project_name = project_name,
	  		Genes = Genes,
	  		Gene_Buffer = Gene_Buffer,
	  		SNPs = SNPs,
	  		SNP_names = SNP_names,
	  		SNP_buffer = SNP_buffer
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
		Array[String]? SNPs
		Array[String]? SNP_names
		Int? SNP_buffer
	}
	String output_file = "~{project_name}.genomic_intervals.bed"
	Array[String] SNPs_to_use = select_first([SNPs, []])
	Array[String] SNP_names_to_use = select_first([SNP_names, []])
	Int SNP_buffer_to_use = select_first([SNP_buffer, 0])

	command <<<
	set -euxo pipefail

	gunzip -c ~{MANE_GENCODE} > MANE.GRCh38.gtf

	python3 <<CODE
	import pandas as pd

	f.open("~{output_file}","w")

	mane_df = pd.read_csv("MANE.GRCh38.gtf",
		sep='\t',
		index_col=False,
		comment='#',
		header=None,
		names = ["chrom","source","feature","start","end","score","strand","frame","attributes"])

	genes = "~{sep=' ' Genes}".split()
	
	# Loop through genes
	for gene in genes:
		gene_df = mane_df[mane_df["attributes"].str.contains(f"gene_name \"{gene}\"", na=False)]
		gene_df = gene_df[gene_df['feature'] == "gene"]
		f.write((gene_df['chrom']) + "\t" + str(int(gene_df['start']) - ~{Gene_Buffer}) + "\t" + str(int(gene_df['end']) + ~{Gene_Buffer}) + "\t" + gene)


	CODE
	# # Loop through all of the genes in a list to find relevant genomic ranges
	# for gene in ~{sep=' ' Genes}; do
	# 	awk -v gene="$gene" '
	# 		$3 == "gene" {
	# 			n = split($9, attrs, ";")
	# 			for (i = 1; i <= n; i++) {
	# 				if (attrs[i] ~ /^[[:space:]]*gene_name[[:space:]]/) {
	# 					gsub(/^[[:space:]]*gene_name[[:space:]]*"/, "", attrs[i])
	# 					gsub(/"[[:space:]]*$/, "", attrs[i])
	# 					if (attrs[i] == gene) {
	# 						print $1 "\t" ($4 - 1) "\t" $5 "\t" gene
	# 					}
	# 				}
	# 			}
	# 		}
	# ' MANE.GRCh38.gtf >> ~{output_file}
	# done

	# # Loop through all SNPs in a list
	# if [[ ~{length(SNPs_to_use)} -gt 0 ]]; then

	# 	i=0
	# 	for snp in ~{sep=' ' SNPs_to_use}; do
	# 		chrom=$(echo "$snp" | cut -d: -f1)
	# 		pos=$(echo "$snp" | cut -d: -f2)

	# 		min=$((pos - SNP_buffer_to_use))
	# 		max=$((pos + SNP_buffer_to_use))

	# 		name=$(echo "~{sep=' ' SNP_names_to_use}" | cut -d' ' -f$((i + 1)))

	# 		echo -e "$chrom\t$min\t$max\t$name" >> ~{output_file}

	# 		i=$((i + 1))
	# 	done
	# fi
	>>>
	runtime {
		docker: "vanallenlab/g2c_pipeline"
	}
	output {
		File out1 = "~{output_file}"
	}
}