import os
def find_variant(file):
    vf=open(file).readlines()
    if not vf:
        return {}
    vardict={}
    pos=[]
    type=[]
    rep=[]
    linter=[]
    init=False
    inter=[]
    for line in vf:
        lsplit=line.split('\t')
        if ':' in lsplit[0] :
            if init:
                linter=correct_inter(inter) 
                vardict[nom]=[lenght,pos,type,linter,rep]
                pos=[]
                type=[]
                rep=[]
                linter=[]
                inter=[]
            dsplit=lsplit[0].split('_')
            nom=dsplit[0]
            lenght=dsplit[1][:-2]
            init=True
        else:
            pos.append(lsplit[1])
            type.append(lsplit[0])
            
            if lsplit[0]=='D':
                inter.append(1)
                rep.append(lsplit[2][2])
            elif lsplit[0]=='I':
                inter.append(-1)
                rep.append('-')
            else:
                inter.append(0)
                rep.append(lsplit[2][1])
            
    linter=correct_inter(inter)        
    vardict[nom]=[lenght,pos,type,linter,rep]
    return vardict


def correct_inter(inter):
    inter.reverse()
    count=0
    lst=[]
    for i in inter:
        count+=i
        lst.append(count)
    lst.reverse()
    return lst


def cor_sam(file,varfile,ffai,outfile):
    header = open(ffai, 'r').readlines()[0].split('\t')[0]
    sam = open(file, 'r').readlines()
    dic=find_variant(varfile)
    n_sam=open(outfile,'w')
    if extract_spefile(file)[0:2]=='WT':
        for l in sam:
            n_sam.write(l.replace('ref',header))
        return
    for block in range(int(len(sam)/4)):
        if sam[4*block][0]=='a':
            s1,s2=correct_ref(sam[4*block+1],sam[4*block+2],dic[extract_spefile(file)],header)
            n_sam.write('a\n')
            n_sam.write(s1+'\n')
            n_sam.write(s2+'\n\n')
    n_sam.close()
    #os.remove(file)   


def correct_ref(s1,s2,dic,header):
    tab=list(filter(None,s1.split(' ')))
    tab2=list(filter(None,s2.split(' ')))
    space=' '
    if tab2[1]=='S1_6904':
        pass
    lst,stp_pos,first_inter=find_inter(int(tab[2]),int(tab[3]),dic) 
    rep =tab[6]
    rep2=tab2[6]
    if len(lst)>0:
        for pos in lst:
            rep,rep2,cp=correct_seq(rep,int(pos[0]),stp_pos,int(pos[2]),pos[1],pos[3],rep2,first_inter)
        if cp:
            stp_pos+=1

        l_frag=int(tab[3])+(next_inter(lst[0])-lst[-1][2])
        if l_frag!=len(rep.replace('-','').replace('\n','')):
            print(l_frag,len(rep.replace('-','').replace('\n','')),tab[3],'\n',tab[2],lst)
            raise  NameError('s1 error in length') 
        if int(tab2[3])!=len(rep2.replace('\n','').replace('-','')):
            print(tab2[3],len(rep2.replace('\n','').replace('-','')),'\n',tab2[0],tab2[1],tab2[2],tab2[3],tab2[4])
            raise  NameError('s2 error in lenght') 
        s1=('{0}{7}{1}{7}{2}{7}{3}{7}{4}{7}{5}{7}{6}'.format(tab[0],header,str(stp_pos),str(l_frag),tab[4],str(dic[0]),rep,space)).replace('\n','')
        s2=('{0}{7}{1}{7}{2}{7}{3}{7}{4}{7}{5}{7}{6}'.format(tab2[0],tab2[1],tab2[2],tab2[3],tab2[4],tab2[5],rep2,space)).replace('\n','')
    else:
        s2=s2.replace('\n','')

        s1=('{0}{7}{1}{7}{2}{7}{3}{7}{4}{7}{5}{7}{6}'.format(tab[0],header,str(stp_pos),tab[3],tab[4],str(dic[0]),tab[6],space)).replace('\n','')
    return s1,s2
        
    
def prev_inter(tup):
    rep=int(tup[2])
    if tup[1]=='D':
        rep-=1
    elif tup[1]=='I':
        rep+=1
    return rep

def next_inter(tup):
    rep=int(tup[2])
    if tup[1]=='D':
        rep+=1
    elif tup[1]=='I':
        rep-=1
    return rep

def find_inter(a,b,dic):
    lst=[]
    i=0
    stp_pos=10000000000000000
    first_inter=0
    if a+prev_inter((0,dic[2][(0)],dic[3][(0)]))<int(dic[1][0]):
        while(i <= len(dic[1])-1 and stp_pos==10000000000000000):
            if int(dic[1][-(i+1)])>a+prev_inter((0,dic[2][-(i+1)],dic[3][-(i+1)])):
                stp_pos=a+prev_inter((0,dic[2][-(i+1)],dic[3][-(i+1)]))
                elem=(0,dic[2][-(i+1)],dic[3][-(i+1)])
                first_inter=prev_inter(elem)
            i+=1
        if stp_pos==10000000000000000:
            raise('error in start_position')
    else:
        stp_pos=a+dic[3][(0)]     
    for idx,p in enumerate(dic[1]):
        if int(p)>stp_pos and int(p)<stp_pos+b+(prev_inter((0,dic[2][idx],dic[3][idx]))-first_inter):
            try:
                lst.append((int(p),dic[2][idx],prev_inter((0,dic[2][idx],dic[3][idx])),dic[4][idx]))
            except:
                print(idx,dic)
                raise 'ok'
    return lst,stp_pos,first_inter

def correct_seq(seq,pos,st_pos,inter,tup,rep,s2,first_inter):
    res=''
    count=st_pos
    blocks=seq.replace('\n','').split('-')
    make=True
    cp=False
    cp2=False
    for idx,block in enumerate(blocks):
        count+=len(block)
        if count+(inter-first_inter)>pos and make:
            if pos==st_pos+1 and tup=='I':
                cp=True
            count-=len(block)
            posi=(pos)-(count+(inter-first_inter))
            if tup=='I' :
                if s2[len(res)+posi]!='-':
                    if posi!=len(block):
                        res+=block[0:posi]+rep+block[posi+1:]+'-'
                    else:
                        if s2[len(res)+posi-1]!='-':
                            res+=block[:-1]+rep+'-'
                        else:
                            s2=s2[0:len(res)+posi-1]+s2[len(res)+posi:]
                            res+=block[:-1]+'-'
                        
                else:
                    if posi!=len(block):
                        s2=s2[0:len(res)+posi]+s2[len(res)+posi+1:]
                        res+=block[0:posi]+block[posi+1:]+'-'
                        
                    else:
                        s2=s2[0:len(res)+posi]+s2[len(res)+posi+1:]
                        res+=block[:-1]+'-'
            elif tup=='D':
                if posi!=len(block):
                    s2=correct_s2(s2,len(res)+posi)
                    res+=block[0:posi]+rep+block[posi:]+'-'
                else:
                    s2=correct_s2(s2,len(res)+posi)
                    res+=block+rep+'-'
            elif tup=='S':
                if posi!=len(block):
                    res+=block[0:posi]+rep+block[posi+1:]+'-'
                else:
                    res+=block[:-1]+rep+'-'
            count+=len(block)
            make=False  
        else:
            res+=block+'-'
    if make:
        raise NameError('one mutation did not be errased')
    if len(res[:-1])!=len(s2.replace('\n','')):
        raise  NameError('s1 s2 different length',len(res[:-1]),len(s2.replace('\n','')))
    if cp:
        if seq[0]!='-' and res[0]=='-':          
            cp2=True
        elif seq[0]=='-'and res[0:2]=='--':
            cp2=True
    return res[:-1],s2,cp2        

def extract_spefile(file):
    temp = file.split('pbsim2')[1].split('_')[1:3]
    return temp[1]+temp[0]
    
def correct_s2(s2,posi):
    return s2[0:posi]+'-'+s2[posi:]


def correct_doubleinsert(a,b):
    i=0
    while(i < len(a)):
        if a[i]==b[i]==['-']:
            a=a[0:i]+a[i+1:]+'-'
            b=b[0:i]+b[i+1:]+'-'
            i-=1
        i+=1
    return a,b





def count_pos():
    pass

def do_something(filein, fileout):
    cor_sam(filein[1]+'_0001.maf',filein[2],filein[3],fileout[0])
   



do_something(snakemake.input, snakemake.output) 
#cor_sam('Workflow/result/pbsim2/simu_covid_V0_150_0.99_0001.maf','Workflow/data/covid_variant_file.txt','Workflow/data/ref_covid.fasta.fai','Workflow/result/pbsim2/simu_covid_V0_150_0.99_0001_corrcted.maf')
