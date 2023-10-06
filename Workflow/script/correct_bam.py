import pysam
def extract_reads(tools, minimap):
    """
    Extract all the reads from several BAM files.

    Parameters
    ----------
    *args : str
        x BAM filename
    perfect : str
        Name of the BAM file containing the perfect alignment
        
    Return
    ------
    - 1/ Dictionnary storing all reads contained in the "perfect" BAM file.
        --> Keys : "Number from the name of the read"
        --> Values : Read corresponding to the number (pysam.AlignedSegment object)
    
    - 2/ Dictionnary storing all reads for each BAM file from the different mapping tools
        --> Keys : file1 , file2, ..., filen
        --> Values : Another dictionnary containing all reads from one file in the same format that the dictionnary create for the "perfect" file.
    
    - 3/ Length of the sequence reference (int)
    """
    reads_perfect = {}
    filemini =open(minimap, "r").readlines()
    for r in filemini:
        opt=r.split('\t')
        if r[0]!='@':
            if opt[1]!=4:
                reads_perfect[opt[0]] = opt

    all_files = {}
    c = 1
    for tool in tools:
        toomfile =open(tool, "r").readlines()
        newsam=open(tool.replace('.sam','_corrected.sam'),'w')
        for r in toomfile:
            if r[0]=='@':
                newsam.write(r)
            else:
                opt=r.split('\t')
                if int(opt[1])>2000:
                    continue
                if opt[1]!='4':
                    if opt[0] in reads_perfect.keys() : 
                        min_r=reads_perfect[opt[0]]
                        if len(min_r[6])!=len(opt[6]):
                            print(":()")
                        begin_r=int(opt[3])
                        cig=opt[5]
                        cig_pos=findcig(cig)
                        m_cig=findcig(min_r[5])
                        m_begin=int(min_r[3])
                        npos=begin_r
                        if m_begin-m_cig[0]!=begin_r-cig_pos[0]:
                            print(":'()")
                        if m_cig[0]>0 or m_cig[1]>0:
                            cig=cigttuple(cig)
                            if int(m_cig[0])>int(cig_pos[0]):
                                cig,npos=new_cig(cig,m_cig[0])
                                npos=begin_r+npos
                            if int(m_cig[1])>int(cig_pos[1]):
                                cig,_=new_cig(cig,(m_cig[1]),True)
                            nbopt=len(opt)-1
                            for i in range(len(opt)):
                                toadd=opt[i]+'\t'
                                if i==3:
                                    toadd=str(npos)+'\t'
                                if i==5:
                                    toadd=''
                                    for lcig in cig:
                                        for elem in lcig:
                                            toadd=toadd+str(elem)
                                    toadd=toadd+'\t'
                                if i==nbopt:
                                    toadd=opt[i]
                                newsam.write(toadd)
                        else:   
                            newsam.write(r)
                    else:       
                        newsam.write(r)
                else:
                    newsam.write(r)
    return 

def findcig(tup):
    rep1=0
    rep2=0
    for a in range(len(tup)):
        if not tup[a].isnumeric():
            if tup[a]=='S':
                rep1=int(tup[:a])
                break
            else:
                break
    if tup[-1]=='S':
        for i,x in enumerate(tup[:-1][::-1]):
            if not x.isnumeric():
                rep2=int(tup[len(tup)-(i+1):len(tup)-1])
                break
    return (rep1,rep2)

def new_cig(cig,correct,end=False):
    cor=correct
    init=False
    #todel=len(str(correct[1]))
    ls=cig.copy()
    res=[]
    npos=0
    if end:
        ls=ls[::-1]
    if ls[0][1]!='S':
        ls.insert(0,[0,'S'])
    count=cor
    for tup in ls:
        if tup[1]=='S' and init==False:
            res.append([cor,'S'])
            count-=tup[0]
            init=True
        else:
            if count>=0:
                if not tup[1]=='I':
                    npos+=tup[0]

                if not tup[1]=='D':
                    if count-tup[0]<0:
                        res.append([tup[0]-count,tup[1]])
                        count=-1
                        npos-=count
                    elif count-tup[0]==0:
                        count=-1
                    else:
                        count=count-tup[0]  
            else:
                res.append(tup)
    if end:
        res.reverse()
    return res,npos

def cigttuple(tup):
    cig=[]
    d=0
    for a in range(len(tup)):
        if tup[a].isnumeric():
            pass
        else:
            cig.append([int(tup[d:a]),tup[a]])
            d=a+1
    return cig



def do_something():
    print(snakemake.input['a'],snakemake.input['b'])
    extract_reads([snakemake.input['a']],snakemake.input['b'])

do_something()
