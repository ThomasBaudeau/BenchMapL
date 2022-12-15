import random

def extract_param_out(outpath):
    rep=outpath.split('/')[2].split('_')
    size=extract_len_specie(rep[0])
    msize=rep[1]
    return ((int(size)/int(msize)))

    

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

def extract_read(file,nbread,comb):
    tableau=open(file,'r').readlines()
    print(nbread)
    print(int((len(tableau)-1)/2))
    lst=random.sample(range(0, int(((len(tableau)-1)/2)-1000)), nbread)
    for i in lst:
        comb+=(tableau[(i*2)])
        comb+=(tableau[(i*2+1)])
    return comb

def do_something(data_path, out_path,param):
    """main function for snakemale

    :param data_path: snakemake.input
    :type data_path: string
    :param out_path: snakemake.output
    :type out_path: string
    """
    comb=open(out_path[0],'a')
    taux=attrib_percent(param['variant']['taux'])
    nb=extract_param_out(out_path[0])
    res=''
    for inp in range(len(data_path)):
        res=extract_read(data_path[inp],int(taux[inp]*(nb/10)),res)
    comb.write(res)
    comb.close()

do_something(snakemake.input,snakemake.output,snakemake.config)


