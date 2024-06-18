def make_bed(a,v):
    ref=open(str(a),'r').readlines()
    bed=open(str(v),'w')
    debut=31
    contig=ref[0].split(' ')[0].split('|')[-1].replace('>','').replace('\n','')
    print(contig)
    seq=''.join(ref[1:]).replace('\n','')
    fin=len(seq)
    if seq[-1]=='A':
        i=-1
        while (seq[i]=='A'):
            i+=-1
        fin=fin+i
    bed.write(contig+'\t'+str(debut)+'\t'+str(fin))


make_bed(snakemake.input,snakemake.output)