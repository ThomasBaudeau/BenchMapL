import functools
import operator
import random

def mutate(nuc):
    nucleotide=['A','T','C','G']
    mut = nucleotide[random.randint(0, 3)]
    while mut==nuc:
        mut = nucleotide[random.randint(0, 3)]
    return mut 

def makerandom(file,wt,input,output,param):
    seq=''.join(wt.readlines()[1:]).replace('\n','')
    name='>'#todo
    nbevent=random.randint(0,2)
    if nbevent==0:   
        seq,name=randominsert(seq, name, pos='random')
    if nbevent == 1:
        seq, name = randomdel(seq, name, pos='random')
    if nbevent == 2:
        seq, name = randomsub(seq, name, pos='random')
    file.write(name)
    file.write(seq)
    file.close()
    return 


def sperandom(file, wt, input, output, param):
    seq = ''.join(wt.readlines()[1:]).replace('\n', '')
    name='>'#todo
    #todo pour extraire les params et faire attention au position des différentes mutations
    if nbevent==0:   
        seq,name=randominsert(seq, name, pos='random')
    if nbevent == 1:
        seq, name = randomdel(seq, name, pos='random')
    if nbevent == 2:
        seq, name = randomsub(seq, name, pos='random')
    file.write(name)
    file.write(seq)
    file.close()
    return 


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
            makerandom(file, wt, input, output, param)
        else:
            sperandom(file, wt, input, output, param)
        
    files1.close()
    files2=open(snakemake.output[0],'a')
    files2.close()
    os.remove(snakemake.input[2]+'_0001.fastq')
    os.remove(snakemake.input[2]+'_0001.maf')
    os.remove(snakemake.input[2]+'_0001.ref')




def randominsert(seq,name,pos='random'):
    if pos == 'random':
        pos=random.randint(0,len(seq)-1)
    mut=mutate(None)
    seq = seq[:pos]+mut+seq[pos:]
    name += str(pos+1)+'I_'
    return seq,name


def randomdel(seq, name, pos='random'):
    if pos == 'random':
        pos = random.randint(0, len(seq)-1)
    seq = seq[:pos-1]+seq[pos:]
    name += str(pos)+'D_'
    return seq, name

def randomsub(seq, name, pos='random'):
    if pos == 'random':
        pos = random.randint(0, len(seq)-1)
    mut=mutate(seq[pos])
    seq = seq[0:pos-1]+mut+seq[pos:]
    name += str(pos)+'S_'
    return seq,name
