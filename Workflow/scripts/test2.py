import pysam



bcf_resu=open('data/variant_file.txt','r')
tab=bcf_resu.readlines()
result={}
posvar=0
for elem in tab:
    data=elem.split('\t')
    result[data[2]]=data[3]
    posvar+=1
bcf_in=pysam.VariantFile("medaka/graphmap_vih_350_0.9_-k#4.vcf")
for rec in bcf_in.fetch():
    if rec.pos in result:
        #incremTP
        posvar-=1
        print(rec.alleles,result[rec.pos])
    else :
        pass #increm FP
for _ in range(posvar):
    pass #increm FN