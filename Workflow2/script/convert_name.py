import re

def openfile(file):
    """open a file 

    :param file: name of a file
    :type file: string
    :return: list of all the lines of the files
    :rtype: list
    """
    return list(open(file).readlines())

def filt_1(x):
    if x[0]=='s':
        return x

def filt_2(x):
    pattern = "^[ATCG]*$"
    if bool(re.match(pattern,x[0:10])) or x[0:2]=='@S':
        return x

def findr(x):
    """filter for find the strand of a reads

    :param x: name of the reads
    :type x: string
    :return: F foward, R reverse
    :rtype: string
    """
    if list(filter(None,x.split(' ')))[4]=='+':
        return 'F'
    else: 
        return 'R'

def do_something():#(data_path, out_path, myparam):
    """main function for rename all the reads in a fasta file
    """
    tab1=list(filter(filt_2,openfile(snakemake.input[0]+'_0001.fastq')))
    tab2=list(filter(filt_1,openfile(snakemake.input[0]+'_0001.maf')))
    files1=open(snakemake.output[0],'a')
    for nb in range(int(len(tab2)/2)):
        temp=list(filter(None,tab2[nb*2].split(' ')))
        files1.write(tab1[2*nb].replace('_','N').replace('\n','')+'_'+temp[2]+'_aligned_0_'+findr(tab2[(nb*2)+1])+'_0_'+temp[3]+'\n')
        files1.write(tab1[(2*nb)+1])
    files1.close()
do_something()#(snakemake.input[0], snakemake.output[0], snakemake.config["param"])