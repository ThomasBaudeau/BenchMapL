import re
from tqdm import tqdm

def openfile(file):
    """open a file 

    :param file: name of a file
    :type file: string
    :return: list of all the lines of the files
    :rtype: list
    """
    with open(file,'r') as myfile:
        return list(myfile.readlines())

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
    print('#############  begin correct name ############# ')
    tab3=openfile(snakemake.input[0]+'_0001.fastq')
    print('Finish Step 1 : fastq open ')
    tab2=list(filter(filt_1,openfile(snakemake.input[0]+'_0001.maf')))
    print('Finish Step 2 : maf filtered ')
    files1=open(snakemake.output[0],'a')
    for nb in tqdm(range(int(len(tab2)/2))):
        temp=list(filter(None,tab2[nb*2].split(' ')))
        files1.write('@'+snakemake.output[0].split('_')[2]+'_'+(tab3[4*nb].replace('_','N').replace('\n','')+'_'+temp[2]+'_aligned_0_'+findr(tab2[(nb*2)+1])+'_0_'+temp[3]+'\n')[1:])
        files1.write(tab3[(4*nb)+1])
        files1.write("+\n")
        files1.write(tab3[(4*nb)+3])
    print("############# end correct name ############# ")
    files1.close()
do_something()#(snakemake.input[0], snakemake.output[0], snakemake.config["param"])