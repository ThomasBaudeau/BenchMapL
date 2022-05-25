import random
import os

def select_reads():
    seq=[]
    random.seed(snakemake.config['seed'])
    files0=open(snakemake.input[0])
    lignes=list(files0.readlines())
    b=int(((len(lignes))/4)-1)
    for i in range(1000):
        a=random.randint(0,b)*4
        seq.append(lignes[a+1].replace('\n',''))
    files0.close()
    return seq


def do_something(data_path, out_path):
    seq=select_reads()
    print(snakemake.input)
    files1=open(snakemake.input[1]+'_aligned_reads.fasta','a')
    for i in range(len(seq)):
        files1.write('\n>ERR3278963_-1000_human_-1_R_-1_-1_-1'+str(i))
        files1.write('\n'+seq[i][0:int(snakemake.input[1].split('_')[2])])
    files1.close()
    files2=open(snakemake.output[0],'a')
    files2.close()
    os.remove(snakemake.input[1]+'_aligned_error_profile')
    os.remove(snakemake.input[1]+'_unaligned_reads.fasta')

do_something(snakemake.input[0], snakemake.output[0])

