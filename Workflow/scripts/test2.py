import pysam

bcf_in=pysam.VariantFile("medaka/graphmap_vih_350_0.9_-k#4.vcf")
for rec in bcf_in.fetch():
    print(rec.alleles,rec.pos)