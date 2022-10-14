import functools
import operator

def makerandom(file,wt,input,output,param):
    seq=''.join(wt.readlines()[1:]).replace('\n','')
    randominsert(seq,pos='random')
    randomdel(seq,pos='random')
    randomsub(seq,pos='random')

def sperandom(file,input,output,param):
    pass


def do_something(data_path, out_path,param):
    """main function for snakemale

    :param data_path: snakemake.input
    :type data_path: string
    :param out_path: snakemake.output
    :type out_path: string
    """
    nbvariant=param['variant']['number']
    spevariant=param['variant']['specify']
    variant_file=open('variant_file.txt','a')
    for variant in range(param['variant']['number']):
        wildtype=open(snakemake.input[1],'r')
        file2=open(snakemake.input[1],'a')
        if len(spevariant)<variant:
            makerandom(file2,wildtype)
        else:
            sperandom
        
    files1.close()
    files2=open(snakemake.output[0],'a')
    files2.close()
    os.remove(snakemake.input[2]+'_0001.fastq')
    os.remove(snakemake.input[2]+'_0001.maf')
    os.remove(snakemake.input[2]+'_0001.ref')


ok=open('Workflow2/data/ref_vih.fasta','r')
print()

randominsert(seq,pos='random')
randomdel(seq,pos='random')
randomsub(seq,pos='random')