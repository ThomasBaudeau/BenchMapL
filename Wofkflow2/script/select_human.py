import random

def select_reads():
    seq=[]
    random.seed(snakemake.config['seed'])
    files0=open(snakemake.input[0])
    lignes=list(files0.readlines())
    b=((len(lignes))/4)-1
    print(b)
    for i in range(1000):
        a=random.randint(0,b)*4
        print(lignes[a])
        seq.append(lignes[a+1].replace('\n',''))
    files0.close()
    return seq



def do_something(data_path, out_path):
    seq=select_reads()
    print(snakemake.input)
    files1=open(snakemake.input[1]+'aligned_reads.fasta','a')
    for i in range(len(seq)):
        files1.write('\n>ERR3278963_human_'+str(i))
        files1.write('\n'+seq[i][0:int(snakemake.input[1].split('_')[2])])
    files1.close()
    files2=open(snakemake.output[0],'a')
    files2.close()

do_something(snakemake.input[0], snakemake.output[0])

