import functools
import operator
import random

def mutate(nuc):
    nucleotide=['A','T','C','G']
    mut = nucleotide[random.randint(0, 3)]
    while mut==nuc:
        mut = nucleotide[random.randint(0, 3)]
    return mut 

def makerandom(file,wt,datapath,variant,tx,fvar):
    seq=''.join(wt[1:]).replace('\n','')
    species='>'+find_species(datapath)+'var@'+str(variant)+'\n'
    name=''
    for _ in range(calcnbmut(len(seq),(tx*2500))):
        nbevent=random.randint(0,2)
        fvar.write(find_species(datapath)+'var@'+str(variant)+'\t')
        if nbevent==0:   
            seq,name=randominsert(seq, name, pos='random')
        if nbevent == 1:
            seq, name = randomdel(seq, name, pos='random')
        if nbevent == 2:
            seq, name = randomsub(seq, name, pos='random')
        fvar.write(name)
    file.write(species)
    file.write(seq)
    file.close()
    return 

def calcnbmut(ln,nb):
    return int(ln/nb)

def randominsert(seq,name,pos='random'):
    if pos == 'random':
        pos=random.randint(0,len(seq)-1)
    mut=mutate(None)
    res = seq[0:pos]+mut+seq[pos:]
    name ='I\t{0}\t({1},{1}{2})\n'
    mut1=seq[pos-1]
    if pos-1<0:
        mut1='.'
    return res,name.format(str(pos+1),mut1,mut)


def randomdel(seq, name, pos='random'):
    if pos == 'random':
        pos = random.randint(0, len(seq)-1)
    nuc=seq[pos-2:pos]
    res = seq[0:pos]+seq[pos+1:]
    name = "D\t{0}\t({1},{2})\n"
    return res, name.format(str(pos),nuc,seq[pos-2])

def randomsub(seq, name, pos='random'):
    if pos == 'random':
        pos = random.randint(0, len(seq)-1)
    mut=mutate(seq[pos])
    res = seq[0:pos]+mut+seq[pos+1:]
    name = "S\t{0}\t({1},{2})\n"
    return res,name.format(str(pos+1),seq[pos],mut)


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
    spevariant=param['variant']['specify']
    resuvariant=open('variant_file.txt','w')
    for idx,variant in enumerate(out_path):
        file2=open(variant,'a')
        if idx==0:
            file2.write(wildtype[0])
            file2.write(''.join(wildtype[1:]).replace('\n',''))
        elif isinstance(spevariant[idx-1],(int,float)):
            makerandom(file2, wildtype, data_path[0], idx,spevariant[idx-1],resuvariant)
        else:
            pass
    file2.close()
    resuvariant.close()


do_something(snakemake.input,snakemake.output,snakemake.config)
