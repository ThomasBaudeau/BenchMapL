import random
import os
import time

def openfile(file):
    """open a file 

    :param file: name of a file
    :type file: string
    :return: list of all the lines of the files
    :rtype: list
    """
    with open(file,'r') as myfile:
        return list(myfile.readlines())


def select_reads(nb):
    """select reads and add them in a list

    :return: list of reads
    :rtype:list
    """
    seq=[]
    score=[]
    random.seed(snakemake.config['seed'])
    files0=open(snakemake.input[0])
    lignes=list(files0.readlines())
    b=((len(lignes))/4)-1
    for _ in range(round((snakemake.config["contamination"]/100)*nb)):
        a=random.randint(0,b)*4
        seq.append(lignes[a+1].replace('\n',''))
        score.append(lignes[a+3].replace('\n',''))
        print(':()')
    files0.close()
    print(len(seq))

    return seq,score



def do_something(data_path, out_path):
    """main function for snakemale

    :param data_path: snakemake.input
    :type data_path: string
    :param out_path: snakemake.output
    :type out_path: string
    """
    random.seed(int(snakemake.config['seed']))
    seq,score=select_reads(len(list(openfile(snakemake.input[1])))/4)
    files1=open(snakemake.input[1],'a')
    val=3
    if snakemake.input[1].split('_')[3]=='variantreads1.fastq':
        val=1
    for i in range(len(seq)):
        print(snakemake.input[1].split('_')[val])
        files1.write('@ERR3278963_-1000_human_-1_R_-1_-1_-1'+str(i))
        files1.write('\n'+seq[i][0:int(snakemake.input[1].split('_')[val])])
        files1.write('\n+\n'+score[i][0:int(snakemake.input[1].split('_')[val])]+'\n')
    files1.close()
    os.rename(snakemake.input[1],out_path)


do_something(snakemake.input[0], snakemake.output[0])

