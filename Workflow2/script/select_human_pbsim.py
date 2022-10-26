import random
import os
def select_reads():
    """select reads and add them in a list

    :return: list of reads
    :rtype:list
    """
    seq=[]
    random.seed(snakemake.config['seed'])
    files0=open(snakemake.input[0])
    lignes=list(files0.readlines())
    b=((len(lignes))/4)-1
    for _ in range(1000):
        a=random.randint(0,b)*4
        seq.append(lignes[a+1].replace('\n',''))
    files0.close()
    return seq



def do_something(data_path, out_path):
    """main function for snakemale

    :param data_path: snakemake.input
    :type data_path: string
    :param out_path: snakemake.output
    :type out_path: string
    """
    seq=select_reads()
    files1=open(snakemake.input[1],'a')
    for i in range(len(seq)):
        files1.write('\n>ERR3278963_-1000_human_-1_R_-1_-1_-1'+str(i))
        files1.write('\n'+seq[i][0:int(snakemake.input[1].split('_')[3])])
    files1.close()
    files2=open(snakemake.output[0],'a')
    files2.close()
    os.remove(snakemake.input[2]+'_0001.fastq')
    os.remove(snakemake.input[2]+'_0001.maf')
    os.remove(snakemake.input[2]+'_0001.ref')
do_something(snakemake.input[0], snakemake.output[0])

