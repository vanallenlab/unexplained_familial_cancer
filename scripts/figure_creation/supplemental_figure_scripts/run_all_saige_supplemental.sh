python3 saige_gene_supplemental.py --cancer breast --criteria "T1;T2" --max_AF 001 --genes BRAT1 --color "#FF6F61" --out breast.qq.pdf
python3 saige_gene_supplemental.py --cancer thyroid --criteria "T1;T2;T3;T4" --max_AF 0001 --genes TSTD2 --color "#9EDAE5" --out null
python3 saige_gene_supplemental.py --cancer kidney --criteria "T1;T2" --max_AF 001 --genes ZNF346 --color "#98DF8A" --out null
python3 saige_gene_supplemental.py --cancer kidney --criteria "T1;T2;T3" --max_AF 001 --genes ZNF346,OLFML2B --color "#98DF8A" --out null
python3 saige_gene_supplemental.py --cancer breast_patient_and_family_prs --criteria "T1" --max_AF 0001 --genes CCDC27 --color "#FF6F61" --out null
python3 saige_gene_supplemental.py --cancer bladder --criteria "T1" --max_AF 0001 --genes MYO1A --color "#AEC7E8" --out null
