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
    lmut=[]
    #tirer au sort toute les mutations puis les appliquer dans l'ordre decroissant
    fvar.write('V'+str(variant-1)+find_species(datapath)+'_'+str(len(seq))+':\n')
    for pos in random.sample(range(100, len(seq)-100), calcnbmut(len(seq),tx)):
        
        nbevent=random.randint(0,2)
        if nbevent==0:   
            name=randominsert(seq, pos)
        if nbevent == 1:
            name = randomdel(seq, pos)
        if nbevent == 2:
            name = randomsub(seq, pos)
        lmut.append(name)
    lmut=correct_mut(sorted(lmut,key=lambda ok : ok[1],reverse=True))
    makemute(file,seq,species,lmut,fvar)

    return 


def correct_mut(lmut):
    for x in range(len(lmut)-2):
        a=lmut[x]
        b=lmut[x+1]
        if int(a[1])-1==int(b[1]):
            if a[0]=='I' and b[0]=='D':
                ori=a[2][1]
                lmut[x]=('S',a[1],'({0},{1})'.format(ori,mutate(ori)))
    return lmut


def calcnbmut(ln,nb):
    return int(ln*(nb/100))

def makemute(file,seq,species,lmut,fvar):
    file.write(species)
    for tup in lmut:
        pos=tup[1]
        if tup[0]=='I':
            seq=seq[0:pos]+tup[2][4]+seq[pos:]
        if tup[0]=='D':
            seq=seq[0:pos]+seq[pos+1:]
        if tup[0]=='S':
            seq=seq[0:pos]+tup[2][3]+seq[pos+1:]
        fvar.write(tup[0]+'\t'+str(pos)+'\t'+tup[2]+'\n')
    file.write(seq)
    file.close()



def randominsert(seq,pos='random'):
    if pos == 'random':
        pos=random.randint(0,len(seq)-1)
    mut=mutate(None)
    mut1=seq[pos-1]
    if pos-1<0:
        mut1='.'
    return ('I',pos,'({0},{0}{1})'.format(mut1,mut))


def randomdel(seq, pos='random'):
    if pos == 'random':
        pos = random.randint(0, len(seq)-1)
    nuc=seq[pos-1:pos+1]
    #todo position 0 voir au dessus
    return ("D",pos,"({0},{1})".format(nuc,seq[pos-1]))

def randomsub(seq, pos='random'):
    if pos == 'random':
        pos = random.randint(0, len(seq)-1)
    mut=mutate(seq[pos])
    return ("S",pos,"({0},{1})".format(seq[pos],mut))


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
    random.seed(int(param['seed']))
    wildtype=open(data_path[0],'r').readlines()
    spevariant=param['variant']['specify']
    resuvariant=open(data_path[0].split('_')[1].replace('.fasta','_')+'variant_file.txt','w')
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
