import hail as hl

hl.init(
    gcs_requester_pays_configuration="cordial-diamond-9893"
)
hl.default_reference("GRCh38")

mt = hl.read_matrix_table(
    "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/exome/splitMT/hail.mt"
)

intervals = []

with open("comprehensive_genes.coordinates.bed") as f:
    for line in f:
        chromosome, start, end, gene = line.rstrip("\n").split("\t")
        intervals.append(f"{chromosome}:{start}-{end}")

print(f"Extracting {len(intervals)} intervals", flush=True)

mt_filtered = hl.filter_intervals(
    mt,
    [hl.parse_locus_interval(x) for x in intervals]
)

hl.export_vcf(
    mt_filtered,
    "gs://dataproc-staging-wb-cordial-diamond-9893/comprehensive_cpgs.vcf.bgz"
)