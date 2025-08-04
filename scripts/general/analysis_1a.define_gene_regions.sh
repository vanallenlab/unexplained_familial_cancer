gunzip -c /Users/noah/Downloads/gencode.v48.annotation.gtf.gz | awk '$3 == "gene"' | grep -E 'lncRNA|protein_coding_gene' | cut -f1,4,5 > part1.bed

gunzip -c /Users/noah/Downloads/gencode.v48.annotation.gtf.gz \
  | awk '$3 == "gene" && $0 ~ /lncRNA|protein_coding/ {
      n = split($0, fields, ";");
      for (i = 1; i <= n; i++) {
          if (fields[i] ~ /gene_name/) {
              gsub(/.*gene_name "|"/, "", fields[i]);
              print fields[i];
              break;
          }
      }
  }' > part2.bed

  paste part1.bed part2.bed > analysis_1a.bed
  rm part1.bed part2.bed