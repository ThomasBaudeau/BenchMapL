import random

def extract_param_out(outpath,param):
    rep=outpath.split('/')[2].split('_')
    size=extract_len_specie(rep[0])
    msize=rep[1]
    return ((int(size)/int(msize))*param['number'])

    

def extract_len_specie(species):
    rep=open('data/ref_'+species+'.fasta','r')
    return len(''.join(rep.readlines()[1:]).replace('\n',''))


def attrib_percent(list):
    tot=0
    for i in list:
        tot+=i
    if tot>100:
        raise('wrong input in config2 (taux)')
    else:
        list.insert(0,100-tot)
    return list

def extract_perfect(rname):
    return int(rname.split('_')[1].replace('S1N',''))
def extract_name(rname):
    return rname.split('_')[0].replace('@','')

def extract_read(file,nbread,comb,pf):
    tableau=open(file,'r').readlines()
    print(tableau)
    perfect_file=file.replace('result/pbsim2','perfect_sam').replace('_reads.fastq','.sam')
    file2=open(perfect_file,'r').readlines()
    print(len(tableau),((len(tableau)-1)/4),nbread)
    while (((len(tableau)-1)/4)<nbread):
        nbread-=5
    lst=random.sample(range(0, int(((len(tableau)-1)/4))), nbread)
    for i in lst:
        
        name=tableau[(i*4)]
        pf+=file2[extract_perfect(name)-1].replace('S1_',extract_name(name)+'_S1N')
        comb+=(name)
        comb+=(tableau[(i*4+1)])
        comb+=(tableau[(i*4+2)])
        comb+=(tableau[(i*4+3)])
    return comb,pf

def do_something(data_path, out_path,param):
    """main function for snakemale

    :param data_path: snakemake.input
    :type data_path: string
    :param out_path: snakemake.output
    :type out_path: string
    """
    comb=open(out_path[0],'a')
    combpf=open(out_path[1],'a')
    taux=attrib_percent(param['variant']['taux'])
    nb=extract_param_out(out_path[0],param)
    res=''
    pf=''
    for inp in range(len(data_path)):
        res,pf=extract_read(data_path[inp],int(taux[inp]*(nb/100)),res,pf)
    comb.write(res)
    combpf.write(pf)
    comb.close()
    combpf.close()

do_something(snakemake.input,snakemake.output,snakemake.config)



