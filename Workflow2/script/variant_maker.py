import functools
import operator
import random

def mutate(nuc):
    nucleotide=['A','T','C','G']
    mut = nucleotide[random.randint(0, 3)]
    while mut==nuc:
        mut = nucleotide[random.randint(0, 3)]
    return mut 

def makerandom(file,wt,datapath,variant):
    seq=''.join(wt[1:]).replace('\n','')
    name='>'+find_species(datapath)+'var@'+str(variant)+'@@'
    nbevent=random.randint(0,2)
    if nbevent==0:   
        seq,name=randominsert(seq, name, pos='random')
    if nbevent == 1:
        seq, name = randomdel(seq, name, pos='random')
    if nbevent == 2:
        seq, name = randomsub(seq, name, pos='random')
    name+='xxx\n'
    file.write(name)
    file.write(seq)
    file.close()
    return 


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


def sperandom(file, wt, input_sn, output_sn, param):
    seq = ''.join(wt[1:]).replace('\n', '')
    name='>'+find_species(input_sn)+'var@'+str(param)+'\n'
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

def find_species(input_sn):
    name= input_sn.split('_')[1]
    name=name.split('.')[0]
    return name 


def do_something(data_path, out_path,param):
    """main function for snakemale

    :param data_path: snakemake.input
    :type data_path: string
    :param out_path: snakemake.output
    :type out_path: string
    """
    wildtype=open(data_path[0],'r').readlines()
    nbvariant=param['variant']['number']
    spevariant=param['variant']['specify']
    variant_file=open('variant_file.txt','a')
    for variant in range(len(out_path)):
        file2=open(out_path[variant],'a')
        # if len(spevariant)<variant:
        #     makerandom(file2, wildtype, data_path, out_path, param)
        # else:
        #     sperandom(file2, wildtype, data_path, out_path, param)
        if variant!=0:
            makerandom(file2, wildtype, data_path[0], variant)
        else:
            file2.write(wildtype[0])
            file2.write(''.join(wildtype[1:]).replace('\n',''))
        file2.close()



do_something(snakemake.input,snakemake.output,snakemake.config)
