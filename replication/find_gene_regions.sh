gzip -dc MANE.GRCh38.v1.3.ensembl_genomic.gtf.gz \
| awk -F'\t' -v OFS='\t' -v buffer=2000 '
    BEGIN {
        while ((getline gene < "comprehensive_cpgs.ufc_replication.list") > 0)
            wanted[gene] = 1
    }

    $3 == "gene" {
        gene = $9
        sub(/^.*gene_name "/, "", gene)
        sub(/".*$/, "", gene)

        if (gene in wanted) {
            print $1, $4-buffer, $5+buffer, gene
        }
    }
' \
| sort -k4,4 \
| bgzip -c \
> comprehensive_genes.coordinates.bed.gz
