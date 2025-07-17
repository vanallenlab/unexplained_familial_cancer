
# Subset MANE to genes of interest
zcat MANE.GRCh38.v1.3.ensembl_genomic.gtf.gz \
| fgrep -wf "gencode.v47.autosomal.protein_coding.genes.list" \
| awk '$3 == "gene"' \
> MANE.GRCh38.v1.3.ensembl_genomic.subsetted.gtf

# Extract coordinates from subsetted GTF, pad by ±2kb to be conservative, and use BEDTools to merge these coordinates into unique regions

cat MANE.GRCh38.v1.3.ensembl_genomic.subsetted.gtf \
| awk -F'\t' '
  $3 == "gene" {
    gene_name = $9
    gsub(/.*gene_name "/, "", gene_name)
    gsub(/".*/, "", gene_name)
    print $1, $4, $5, gene_name
  }' \
| sort -Vk1,1 -k2,2n -k3,3n \
> genes.coordinates.bed

