import re
from upsetplot import from_contents
from upsetplot import UpSet
import matplotlib


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam
import tqdm as tq
import os


os.environ['XDG_SESSION_TYPE'] = 'x11'
matplotlib.use('Agg')

class result:
    """class for keept the different result from a bam
    """
    def __init__(self):
       
        self.cor=0
        self.cor_5=0
        self.cor_10=0
        self.cor_20=0
        self.mapped=0
        self.unmapped=0
        self.group1=[]
        self.group2=[]
        self.missaligned=0

    def setname(self,name):
        """set name

        :param name: name of the result
        :type name: str
        """
        self.name=name

    def corpercent(self,cor):
        """percent of cor reads

        :param cor: number of cor reads
        :type cor: int
        :return: percent of correct reads
        :rtype: int
        """
        return 0 if self.mapped==0 else (cor/(self.mapped))*100

    def mappercentU(self):
        """return percent of unmapped reads

        :return: percent of unmapped reads
        :rtype: int
        """
        return ((self.unmapped)/(self.mapped+self.unmapped)*100)

    def mappercentM(self):
        """return percent of mapped reads

        :return: percent of mapped reads
        :rtype: int
        """
        return ((self.mapped)/(self.mapped+self.unmapped)*100)

    def addGroup1(self,add):
        """add reads in group1 reads (unmapped reads)

        :param add: id of the unmapped reads
        :type add: str
        """
        self.group1.append(add)

    def addGroup2(self,add):
        """add reads in group2 reads (misslocalised reads)

        :param add: id of the reads
        :type add: str
        """
        self.group2.append(add)

    def missalign(self):
        """increm missaligned
        """
        self.missaligned+=1

    def increm_cor(self,num=0):
        """increm cor

        :param num: num of the cor, defaults to 0
        :type num: int, optional
        """
        if num not in [0,5,10,20]:
            raise
        elif num==0:
            self.cor+=1
        elif num==5:
            self.cor_5+=1
        elif num==10:
            self.cor_10+=1
        elif num==20:
            self.cor_20+=1

    def get_max_cor(self):
        return self.cor+self.cor_5+self.cor_10+self.cor_20
            

    def increm_mapped(self):
        """increment mapped
        """
        self.mapped+=1

    def increm_unmapped(self):
        """increm mapped
        """
        self.unmapped+=1

class resu_bcf:
    def __init__(self):
       
        self.TP=0
        self.TN=0
        self.FP=0
        self.FN=0
        self.Pr=0
        self.rec=0
        self.f1=0
        self.group=[]
        self.name=''

    def setname(self,name):
        """set name

        :param name: name of the result
        :type name: str
        """
        self.name=name
    
    def addgroup(self,name):
        self.group.append(name)

    def setTP(self):
        self.TP+=1
    def setTN(self):
        self.TN+=1
    def setFP(self):
        self.FP+=1
    def setFN(self):
        self.FN+=1

    def setPr_rec(self):
        if (self.TP+self.FP)==0:
            self.Pr=0
        else:
            self.Pr=self.TP/(self.TP+self.FP)
        if (self.TP+self.FN)==0:
            self.rec=0
        else:
            self.rec=self.TP/(self.TP+self.FN)

    def set_F1(self):
        if self.rec+self.Pr==0:
            self.f1=0
        else:
            self.f1=2*(self.rec*self.Pr)/(self.rec+self.Pr)


def parsenamesimu(text):
    """return expected start position end position and number of reads mapped

    :param text: name of the read
    :type text: str
    :return: expected start position, expected end postion and number of mapped reads
    :rtype: tupple
    """
    try:
        st_pos= int(text.split('_')[2])
        end_pos=int(text.split('_')[7])
    except:
        st_pos='NaN'
        end_pos=''
    return st_pos,end_pos,(st_pos+end_pos)

def read_align(resu,read):
    """increm resu cor for a read

    :param resu: class resu
    :type resu: resu
    :param read: read object from pybam
    :type read: pybam object
    :param threshold: specified cor number, defaults to 0
    :type threshold: int, optional
    :return: result of the operation true if read is increm else false
    :rtype: boolean
    """
    try: 
        st_ref,mbase_ref,end_ref=parsenamesimu(read.query_name)
        st_alg=read.reference_start
        
        end_alg=read.reference_end
        mbase_alg=read.reference_length
        crit1=abs(st_ref-st_alg)
        crit2=abs(end_ref-end_alg)
        for threshold in [0,5,10,20]:
            if (crit1 + crit2)<=threshold :
                resu.increm_cor(threshold)
                return True
        return False 
    except:
        resu.increm_cor(0) 
        return True
    

def countdiff(bamFP):
    """fill the different result for the input bam

    :param bamFP: bam file in pybam format
    :type bamFP: pybam object
    :return: the resu object 
    :rtype: resu
    """
    resu=result()
    for read in bamFP:
        if is_human(read,resu):
            if not read.is_secondary:
                if not read.is_unmapped:
                    resu.increm_mapped()
                    rep=read_align(resu,read)
                    if not rep:
                        resu.addGroup2(read.query_name)
                else:
                    resu.addGroup1(read.query_name)
                    resu.increm_unmapped()
            else:
                pass #todo
    return resu


def find_output(key,name=False,outpath=None):
    """find the output for each groups

    :param key: name of the group
    :type key: str
    :param name: name of the selected output, defaults to False
    :type name: bool or str, optional
    :param outpath: list of all the possible output name, defaults to None
    :type outpath: list, optional
    :return: True if outpout is correct else false
    :rtype: boolean
    """
    for output in outpath:
        if name:
            if output.count(key)!=0 and output.count(name)!=0:
                return output
        else:
            if output.count(key)!=0:
                return True
    return False

def get_label(results):
    """extracts name of all the resu object of the groups

    :param results: list of resu object
    :type results: list
    :return: list of name
    :rtype: list
    """
    labels=[]
    for result in results:
        labels.append(result.name)
    return labels

def is_human(read,resu):
    """tell if read is an human reads and increm missaligned

    :param read: pybam object
    :type read: pybam object
    :param resu: resu object
    :type resu: resu object
    :return: False if missaligned else True
    :rtype: boolean 
    """
    try:
        if read.query_name.split('_')[2]=='human':
            if not read.is_unmapped:
                resu.missalign()
            return False
        return True
    except:
        return True

def get_nbHuman(results):
    """return number of missaligned reads for each result of a groups

    :param results: list of resu object
    :type results: list
    :return: list of the number of missaligned reads for each resu
    :rtype: list
    """
    nb=[]
    for result in results:  
        nb.append(result.missaligned)
    return nb

def get_MapOrNot(results):
    """return number of mapped and unmapped reads for each result of a groups

    :param results: list of resu object
    :type results: list
    :return: lists of the number of mapped and unmapped reads for each resu
    :rtype: list
    """
    data_unmap=[]
    data_map=[]
    for result in results:
        data_unmap.append(result.mappercentU())
        data_map.append(result.mappercentM())
    return data_unmap,data_map

def get_all_cor(results):
    """return list of list of each cor for each result of a groups

    :param results: list of resu object
    :type results: list
    :return: list of list of each cor 
    :rtype: list
    """
    data_cor=[]
    data_cor5=[]
    data_cor10=[]
    data_cor20=[]
    data_corsup=[]
    ttread=[]
    for result in results:
        if result.corpercent(result.get_max_cor())==0:
            data_corsup.append(0)
            data_cor.append(0)
            data_cor5.append(0)
            data_cor10.append(0)
            data_cor20.append(0)
        else:
            data_corsup.append(100-result.corpercent(result.get_max_cor()))
            data_cor.append(result.corpercent(result.cor))
            data_cor5.append(result.corpercent(result.cor_5))
            data_cor10.append(result.corpercent(result.cor_10))
            data_cor20.append(result.corpercent(result.cor_20))
        ttread.append(result.mapped)
    return [data_cor,data_cor5,data_cor10,data_cor20,data_corsup],ttread

def plot_histoMU(results,output):
    """plot histogram of unmapped /mapped reads 

    :param results: list of each resu of the groups
    :type results: list
    :param output: name of the output
    :type output: str
    """
    if output: 
        width = 0.45       # the width of the bars: can also be len(x) sequence
        fig, ax = plt.subplots(figsize=(16, 9))
        plt.rcParams.update({'font.size': 11})
        labels=get_label(results)
        d1,d2=get_MapOrNot(results)
        ax.set_xticklabels(labels, rotation = 45,ha='right')
        ax.bar(labels, d2, width, label='Mapped_reads')
        ax.bar(labels, d1, width, label='Unmapped_reads',bottom=d2)
        ax.set_ylabel('Percentage of reads')
        ax.set_title('Percentage of mapped reads ')
        ax.legend()
        plt.subplots_adjust(left=0.2, bottom=0.2)
        plt.savefig(output,dpi=500,format='pdf')



def plot_histoVar(results,output):
    """plot histogram of unmapped /mapped reads 

    :param results: list of each resu of the groups
    :type results: list
    :param output: name of the output
    :type output: str
    """
    if output:
        labels=[]
        FP=[]
        FN=[]
        TP=[]
        for resu in results:
            labels.append(resu.name)
            if resu.FP<50:
                FP.append(resu.FP)
            else:
                FP.append(40)
            FN.append(resu.FN)
            TP.append(resu.TP)
        x = np.arange(len(labels))  # the label locations
        width = 0.35  # the width of the bars

        fig, ax = plt.subplots()
        rects1 = ax.bar(x - width/3, FP, width, label='FP')
        rects2 = ax.bar(x + width/3, FN, width, label='FN')
        rects3 = ax.bar(x + width/3, TP, width, label='TP')

        # Add some text for labels, title and custom x-axis tick labels, etc.
        ax.set_ylabel('Number')
        ax.set_title('Number of variants by TP FN and FP')
        ax.set_xticks(x, labels)
        ax.legend()

        ax.bar_label(rects1, padding=3)
        ax.bar_label(rects2, padding=3)
        ax.bar_label(rects3, padding=3)

        fig.tight_layout()
        ax.legend()
        plt.subplots_adjust(left=0.2, bottom=0.2)
        plt.savefig(output,dpi=300,format='pdf')



def plot_histoNotHuman(results,output):
    """plot histogram of the number of human reads mapped 

    :param results: list of each resu of the groups
    :type results: list
    :param output: name of the output
    :type output: str
    """
    if output:
        width = 0.45       # the width of the bars: can also be len(x) sequence
        fig, ax = plt.subplots(figsize=(16, 9))
        plt.rcParams.update({'font.size': 11})
        labels=get_label(results)
        ax.set_xticklabels(labels, rotation = 45,ha='right')
        d2=get_nbHuman(results)
        ax.bar(labels, d2, width, label='Human_reads')
        ax.set_ylabel('Percentage of reads')
        ax.set_title('Percentage of mapped contamination reads')
        ax.legend()
        plt.savefig(output,dpi=300,format='pdf')

def common_error_gp1(results,output):
    """plot upsetplot of unmapped reads groups 

    :param results: list of each resu of the groups
    :type results: list
    :param output: name of the output
    :type output: str
    """
    if output:
        plt.clf()
        gp1={}
        for read in results:
            gp1[read.name]=read.group1
        plot1=from_contents(gp1)
        try:
            polt = UpSet(plot1,show_counts=True,element_size=21).plot()
        except:
            pass
        plt.title('Categories of unmapped reads among tools')
        plt.savefig(output,dpi=300,format='pdf')

def common_error_gp2(results,output):
    """plot upsetplot of poorly located reads groups 

    :param results: list of each resu of the groups
    :type results: list
    :param output: name of the output
    :type output: str
    """
    if output:
        plt.clf()
        gp1={}
        try:
            for read in results:
                gp1[read.name]=read.group2
            plot1=from_contents(gp1)
            polt = UpSet(plot1,show_counts=True,element_size=21).plot()
            plt.rcParams.update({'font.size': 11})
            plt.title('Categories of poorly located reads among tools')
            plt.savefig(output)
        except:
            for read in results:
                plt.savefig(output,dpi=300,format='pdf')


def multiple_cor(results,output):
    """plot of the number of reads correctly localised by categories 

    :param results: list of each resu of the groups
    :type results: list
    :param output: name of the output
    :type output: str
    """
    if output:
        plt.clf()
        plt.figure(figsize=(16, 10))
        plt.rcParams.update({'font.size': 16})
        data,ttread = get_all_cor(results)
        columns = get_label(results)
        rows = ['perfect','5b shift','10b shift','20b shift', '>20b shift','# reads']
        values = np.arange(0, 120, 20)
        colors = plt.cm.BuPu(np.linspace(0.15, .65, len(rows)))
        n_rows = len(data)
        index = np.arange(len(columns)) +0.2
        bar_width = 0.5
        y_offset = np.zeros(len(columns))
        y=np.zeros(len(columns))
        cell_text = []
        for row in range(n_rows):
            a=plt.bar(index, data[row], bar_width, bottom=y, color=colors[row])
            y_offset = data[row]
            y=y+data[row]
            #changer la facon de calc le total
            cell_text.append(['%2i' % round(x)  for x in y_offset])
        cell_text.append(['%2i' % round(x)  for x in ttread])
        
        the_table = plt.table(cellText=cell_text,
                            rowLabels=rows, 
                            alpha=1, 
                            rowColours=colors, 
                            colLabels=columns,
                            loc='bottom'
                            )
        the_table.scale(1, 1.7)
        the_table.auto_set_font_size(True)
        #the_table.set_fontsize(15)
        plt.subplots_adjust(left=0.12, bottom=0.19)
        plt.ylabel("Percentage of reads")
        plt.yticks(values, ['%d' % val for val in values])
        # for tx,ypos in zip(ttread,index):
        #     plt.text(ypos,20,tx,color='black')
        plt.xticks([])
        plt.title('Proportion of correctly located read')
        
        plt.savefig(output,dpi=900,format='pdf')

def rsnake(myparam,arg1,arg2=False,arg3=False,arg4=None):
    if arg4!=None:
        return myparam[arg1][arg2][arg3][arg4]
    elif arg2:
        return myparam[arg1][arg2][arg3]
    else:
        return myparam[arg1]


def parsepath(path,myparam):
    """check the command and the tools

    :param path: a path
    :type path: string
    :return: the string name of the command
    :rtype: string
    """
    asplit=path.split('_')
    command=asplit[5][:-4]
    tool=asplit[1][6:]
    for pos in range(len(rsnake(myparam,'param',tool,'command'))):
        
        if rsnake(myparam,'param',tool,'command',pos)==command:

            return rsnake(myparam,'param',tool,'name',pos)
    raise ValueError("Wrong Input")
           

def group_input_param(files,myparam):
    """generate a dict for each group of tools/datasets/command

    :param files: list of bam file from the different tools/dataset
    :type files: list[string]
    :return: dictionnary of the different datasets groups with each files in a list
    :rtype: dict
    """
    group={}
    for specie in rsnake(myparam,"species"):
        for length in rsnake(myparam,'length'):
             for er in rsnake(myparam,'error_rate'):
                for t in rsnake(myparam,'param'):
                    for file in files:
                        if  file.count('_'+specie+'_')==1 and file.count('_'+str(length)+'_')==1  and file.count('_'+str(er)+'_')==1 and file.count(t+'_')==1:
                            try:
                                group[t+'_'+specie+'_'+str(length)+'_'+str(er)+'_'].append(file)
                            except:
                                group[t+'_'+specie+'_'+str(length)+'_'+str(er)+'_']=[]
                                group[t+'_'+specie+'_'+str(length)+'_'+str(er)+'_'].append(file)
    return group

def plot_params(data_path, out_path, myparam):
    """main function

    :param data_path: snakemake.input
    :type data_path: list
    :param out_path: snakemake.output
    :type out_path: list
    :param myparam: config files
    :type myparam: list
    """
    files=list(data_path)
    groups=group_input_param(files,myparam)
    print(groups)
    for key in groups.keys():
        if find_output(key,outpath=out_path): 
            lst=groups[key]
            results=[]
            for tool in lst:
                save = pysam.set_verbosity(0)
                bamFP = pysam.AlignmentFile(tool, "rb")
                pysam.set_verbosity(save)
                name=parsepath(tool,myparam)
                resu=countdiff(bamFP)
                resu.setname(name)
                results.append(resu)
            multiple_cor(results,find_output(key,'bc2', out_path))    
            plot_histoMU(results,find_output(key,'bc1', out_path))
            common_error_gp2(results,find_output(key,'bc3', out_path))

def findname(path):
    """find the name of the tool

    :param path: a path
    :type path: string
    :return: the string between the / and _ character
    :rtype: string
    """
    one=(path.index("/"))
    two=(path.index('_',one))
    return path[one+1:two]           

def findspecies(path):
    """find the name of the tool

    :param path: a path
    :type path: string
    :return: the string between the / and _ character
    :rtype: string
    """
    
    return path.split('_')[1]   



def group_input_simple(files,myparam):
    """generate a dict for each group of tools/datasets

    :param files: list of bam file from the different tools/dataset
    :type files: list[string]
    :return: dictionnary of the different datasets groups with each files in a list
    :rtype: dict
    """
    group={}
    for specie in rsnake(myparam,"species"):
        for length in rsnake(myparam,'length'):
            for er in rsnake(myparam,'error_rate'):
                for file in files:
                    if  file.count('_'+specie+'_')==1 and file.count('_'+str(length)+'_')==1  and file.count('_'+str(er)+'_')==1:
                        try:
                            group[specie+'_'+str(length)+'_'+str(er)+'_'].append(file)
                        except:
                            group[specie+'_'+str(length)+'_'+str(er)+'_']=[]
                            group[specie+'_'+str(length)+'_'+str(er)+'_'].append(file)
    return group

def plot_simple(data_path, out_path, myparam,test=False):
    """main function

    :param data_path: snakemake.input
    :type data_path: list
    :param out_path: snakemake.output
    :type out_path: list
    :param myparam: config files
    :type myparam: list
    """ 
    
    files=list(data_path)
    groups=group_input_simple(files,myparam)
    for key in groups.keys():
        if find_output(key,outpath=out_path): 
            lst=groups[key]
            results=[]
            result_var=[]
            for tool in lst:
                save = pysam.set_verbosity(0)
                bamFP = pysam.AlignmentFile(tool, "rb")
                pysam.set_verbosity(save)
                name=findname(tool)
                resu=countdiff(bamFP)
                file=tool.replace('mapped_reads','medaka').replace('.bam','.vcf')
                resu2=parsevariant(file,parseresuvar(findspecies(file)))
                resu.setname(name)
                results.append(resu)
                result_var.append(resu2)
            common_error_gp1(results,find_output(key,'rl4', out_path))
            common_error_gp2(results,find_output(key,'rl5', out_path))
            multiple_cor(results,find_output(key,'rl2', out_path))    
            plot_histoMU(results,find_output(key,'rl1', out_path))
            plot_histoNotHuman(results,find_output(key,'rl3', out_path))
            upsetplot_var(result_var,find_output(key,'rl6', out_path))
            if test:
                return results

def parsevariant(file,resudic):
    resu=resu_bcf()
    resu.name=findname(file)
    bcf_in=pysam.VariantFile(file)
    posvar=len(resudic.keys())
    for rec in bcf_in.fetch():
        resu.addgroup(str(rec.ref)+'_'+str(rec.pos)+'_'+str(rec.alts))
        if (str((int(rec.pos)-1)) in resudic) or (str((int(rec.pos))) in resudic) or (str((int(rec.pos)-2)) in resudic):
            resu.setTP()
            posvar-=1
        else :
            resu.setFP()
    for _ in range(posvar):
        resu.setFN()
    resu.setPr_rec()
    resu.set_F1()
    return resu



def plotvariant(data_path, out_path, myparam,test=False) :
    files=list(data_path)
    groups=group_input_simple(files,myparam)
    results=[]
    for key in groups.keys():
        if find_output(key,outpath=out_path): 
            lst=groups[key]
            results=[]
            for tool in lst:     
                bcf_in=pysam.VariantFile(tool)
                name=findname(tool)
                resu=parsevariant(bcf_in,parseresuvar())
                resu.setname(findname(tool))
                results.append(resu)
            plot_histoVar(results,out_path[0])




def find_info (file,par):
    param=file[file.index("/")+1:].split('_')
    return {'tname':param[0],'specie':param[1],'lenght':param[2],'er':param[3],'param':param[4],'cov':par['number'],'model':par['model']}


def files_stats(data_path, out_path, myparam):
    """main function

    :param data_path: snakemake.input
    :type data_path: list
    :param out_path: snakemake.output
    :type out_path: list
    :param myparam: config files
    :type myparam: list
    """
    files=list(data_path)
    idx=0
    for file in tq.tqdm(files):
        infos=find_info (file,myparam)
        save = pysam.set_verbosity(0)
        bamFP = pysam.AlignmentFile(file.replace('medaka','mapped_reads').replace('.vcf','.bam'), "rb")
        pysam.set_verbosity(save)
        resu=countdiff(bamFP)
        resu2=parsevariant(file,parseresuvar(findspecies(file)))
        resu3=parsevariant(file.replace('medaka','bcf'),parseresuvar(findspecies(file)))
        infos.update(mapped=resu.mapped,unmapped=resu.unmapped,wrongalign=resu.missaligned,r0=resu.cor,r5=resu.cor_5,r10=resu.cor_10,r20=resu.cor_20,TP=resu2.TP,FN=resu2.FN,FP=resu2.FP,precison=resu2.Pr,recall=resu2.rec,f1score=resu2.f1,bc_TP=resu3.TP,bc_FN=resu3.FN,bc_FP=resu3.FP,bc_precison=resu3.Pr,bc_recall=resu3.rec,bc_f1score=resu3.f1)
        if idx==0:
            df=pd.DataFrame(data=infos,index=[idx])
        else:
            df2=pd.DataFrame(data=infos,index=[idx])
            df=pd.concat([df,df2])
        idx+=1  
    df.to_csv(str(out_path))
    
def upsetplot_var(results,output='test.pdf'):
    if output:
        plt.clf()
        variants={}
        for read in results:
            variants[read.name]=read.group
        plot1=from_contents(variants)
        try:
            polt = UpSet(plot1,show_counts=True,element_size=21).plot()
        except:
            pass
        plt.title('Shared detected variants')
        plt.savefig(output,dpi=300,format='pdf')

def parseresuvar(file=''):
    bcf_resu=open('data/'+file+'_variant_file.txt','r')
    tab=bcf_resu.readlines()
    result={}
    for elem in tab:
        if ':' in elem.replace('\n','')[-1]:
            pass
        else:
            data=elem.split('\t')
            result[data[1]]=data[2]
    return result


#rajouter +1 quand insertion ou délét
