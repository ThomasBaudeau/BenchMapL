import os
def find_variant(file):
    vf=open(file).readlines()
    vardict={}
    pos=[]
    type=[]
    rep=[]
    linter=[]
    init=False
    inter=0
    for line in vf:
        lsplit=line.split('\t')
        if ':' in lsplit[0] :
            if init:
                vardict[nom]=[lenght,pos,type,linter,rep]
            dsplit=lsplit[0].split('_')
            nom=dsplit[0]
            lenght=dsplit[1][:-2]
            init=True
        else:
            pos.append(lsplit[1])
            type.append(lsplit[0])
            
            if lsplit[0]=='D':
                inter-=1
                rep.append(lsplit[2][2])
            elif lsplit[0]=='I':
                inter+=1
                rep.append('-')
            else:
                rep.append(lsplit[2][1])
            linter.append(inter)
            
            
    vardict[nom]=[lenght,pos,type,linter,rep]
    return vardict

def cor_sam(file,varfile,outfile):
    sam=open(file,'r').readlines()
    dic=find_variant(varfile)
    n_sam=open(outfile,'w')
    for block in range(int(len(sam)/4)):
        if sam[4*block][0]=='a':
            s1,s2=correct_ref(sam[4*block+1],sam[4*block+2],dic[extract_spefile(file)])
            n_sam.write('a\n')
            n_sam.write(s1+'\n')
            n_sam.write(s2+'\n\n')
    n_sam.close()
    #os.remove(file)   


def correct_ref(s1,s2,dic):
    tab=list(filter(None,s1.split(' ')))
    tab2=list(filter(None,s2.split(' ')))
    space=' '
    lst=find_inter(int(tab[2]),int(tab[3]),dic)
    if len(lst)>0:
        for pos in lst:
            rep,rep2=correct_seq(tab[6],int(pos[0]),int(tab[2]),int(pos[2]),pos[1],pos[3],tab2[6])
        stp_pos=int(tab[2])-int(lst[-1][2])
        l_frag=int(tab[3])-sum([i[2] for i in lst])
        s1=('{0}{7}{1}{7}{2}{7}{3}{7}{4}{7}{5}{7}{6}'.format(tab[0],tab[1],str(stp_pos),str(l_frag),tab[4],str(dic[0]),rep,space)).replace('\n','')
        s2=('{0}{7}{1}{7}{2}{7}{3}{7}{4}{7}{5}{7}{6}'.format(tab2[0],tab2[1],tab2[2],tab2[3],tab2[4],tab2[5],rep2,space)).replace('\n','')
    else:
        s2=s2.replace('\n','')
        s1=('{0}{7}{1}{7}{2}{7}{3}{7}{4}{7}{5}{7}{6}'.format(tab[0],tab[1],tab[2],tab[3],tab[4],str(dic[0]),tab[6],space)).replace('\n','')
    return s1,s2
        
    



def find_inter(a,b,dic):
    lst=[]
    for idx,p in enumerate(dic[1]):
        if int(p)>a and int(p)<a+b:
            lst.append((int(p),dic[2][idx],dic[3][idx],dic[4][idx]))
    return lst

def correct_seq(seq,pos,st_pos,inter,tup,rep,s2):
    res=''
    count=st_pos
    blocks=seq.split('-')
    make=True
    for idx,block in enumerate(blocks):
        count+=len(block)
        if count>=pos+inter and make:
            count-=len(block)
            posi=pos-count+inter+1
            if tup=='I':
                res+=block[0:posi]+'-'+block[posi+1:]
            if tup=='D':
                s2=correct_s2(s2,len(res)+posi)
                res+=block[0:posi]+rep+block[posi:]
            if tup=='S':
                res+=block[0:posi]+rep+block[posi+1:]
            count+=len(block)
            make=False  
        else:
            res+=block+'-'
    return res[:-1],s2        

def extract_spefile(file):
    return file.split('pbsim2')[1].split('_')[2]
    
def correct_s2(s2,posi):
    return s2[0:posi]+'-'+s2[posi+1:]


def correct_startpos():
    pass

def correct_length():
    pass


def count_pos():
    pass

def do_something(filein, fileout):
    cor_sam(filein[1]+'_0001.maf','variant_file.txt',fileout[0])
   



do_something(snakemake.input, snakemake.output) 

