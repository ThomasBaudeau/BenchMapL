import re
from upsetplot import from_contents
from upsetplot import UpSet
import matplotlib.pyplot as plt
import pysam
import numpy as np


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

    def mappercentU(self,nbreads):
        """return percent of unmapped reads

        :return: percent of unmapped reads
        :rtype: int
        """
        return (((nbreads-1000)-self.mapped)/(nbreads-1000))*100

    def mappercentM(self,nbreads):
        """return percent of mapped reads

        :return: percent of mapped reads
        :rtype: int
        """
        return (self.mapped/(nbreads-1000))*100

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

    def increm_mapped(self):
        """increment mapped
        """
        self.mapped+=1

    def increm_unmapped(self):
        """increm mapped
        """
        self.unmapped+=1


def parsenamesimu(text):
    """return expected start position end position and number of reads mapped

    :param text: name of the read
    :type text: str
    :return: expected start position, expected end postion and number of mapped reads
    :rtype: tupple
    """
    st_pos= int(text.split('_')[1])
    end_pos=int(text.split('_')[6])
    return st_pos,end_pos,(st_pos+end_pos)

def read_align(resu,read,threshold=0):
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
    st_ref,mbase_ref,end_ref=parsenamesimu(read.query_name)
    st_alg=read.reference_start
    end_alg=read.reference_end
    mbase_alg=read.reference_length
    crit1=abs(st_ref-st_alg)
    crit2=abs(end_ref-end_alg)
    if (crit1 + crit2)<=threshold :
        resu.increm_cor(threshold)
        return True
    else:
        return False

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
                    for threshold in [0,5,10,20]:
                        rep=read_align(resu,read,threshold)
                        if rep:
                            break
                        if not rep and threshold==20:
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
    if read.query_name.split('_')[2]=='human':
        if not read.is_unmapped:
            resu.missalign()
        return False
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

def get_MapOrNot(results,nbreads):
    """return number of mapped and unmapped reads for each result of a groups

    :param results: list of resu object
    :type results: list
    :return: lists of the number of mapped and unmapped reads for each resu
    :rtype: list
    """
    data_unmap=[]
    data_map=[]
    for result in results:
        data_unmap.append(result.mappercentU(nbreads))
        data_map.append(result.mappercentM(nbreads))
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
    for result in results:
        data_cor.append(result.corpercent(result.cor))
        data_cor5.append(result.corpercent(result.cor_5))
        data_cor10.append(result.corpercent(result.cor_10))
        data_cor20.append(result.corpercent(result.cor_20))
    return [data_cor,data_cor5,data_cor10,data_cor20]

def plot_histoMU(results,nbreads,output):
    """plot histogram of unmapped /mapped reads 

    :param results: list of each resu of the groups
    :type results: list
    :param output: name of the output
    :type output: str
    """
    if output:
        width = 0.45       # the width of the bars: can also be len(x) sequence
        fig, ax = plt.subplots(figsize=(7, 7))
        plt.rcParams.update({'font.size': 11})
        labels=get_label(results)
        d1,d2=get_MapOrNot(results,nbreads)
        ax.set_xticklabels(labels, rotation = 45)
        ax.bar(labels, d2, width, label='Mapped_reads')
        ax.bar(labels, d1, width, label='Unmapped_reads',bottom=d2)
        ax.set_ylabel('Percentage of reads')
        ax.set_title('Percentage of reads by status and tools')
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
        width = 0.35       # the width of the bars: can also be len(x) sequence
        fig, ax = plt.subplots(figsize=(7, 7))
        plt.rcParams.update({'font.size': 11})
        labels=get_label(results)
        d2=get_nbHuman(results)
        ax.bar(labels, d2, width, label='Human_reads')
        ax.set_ylabel('Percentage of reads')
        ax.set_title('Percentage of non reference\'s reads by tools')
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
        polt = UpSet(plot1,show_counts=True,element_size=21).plot()
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
        plt.figure(figsize=(13, 8))
        plt.rcParams.update({'font.size': 15})
        data = get_all_cor(results)
        columns = get_label(results)
        rows = ['perfect','5pb_shift','10pb_shift','20pb_shift']
        values = np.arange(0, 120, 20)
        colors = plt.cm.BuPu(np.linspace(0.15, .65, len(rows)))
        n_rows = len(data)
        index = np.arange(len(columns)) +0.2
        bar_width = 0.5
        y_offset = np.zeros(len(columns))
        y=np.zeros(len(columns))
        cell_text = []
        for row in range(n_rows):
            plt.bar(index, data[row], bar_width, bottom=y, color=colors[row])
            y_offset = data[row]
            y=y+data[row]
            cell_text.append(['%1i' % (x / 1) for x in y])
        the_table = plt.table(cellText=cell_text,
                            rowLabels=rows, 
                            alpha=1, 
                            rowColours=colors, 
                            colLabels=columns, 
                            loc='bottom')
        the_table.scale(1, 1.8)
        the_table.auto_set_font_size(False)
        the_table.set_fontsize(15)
        plt.subplots_adjust(left=0.12, bottom=0.19)
        plt.ylabel("Percentage of reads")
        plt.yticks(values, ['%d' % val for val in values])
        plt.xticks([])
        plt.title('Proportion of correctly mapped read by tools and error')
        
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
    for key in groups.keys():
        if find_output(key,outpath=out_path): 
            lst=groups[key]
            results=[]
            nbreads=calcnbread(data_path[len(data_path)-1],key)
            for tool in lst:
                save = pysam.set_verbosity(0)
                bamFP = pysam.AlignmentFile(tool, "rb")
                pysam.set_verbosity(save)
                name=parsepath(tool,myparam)
                resu=countdiff(bamFP)
                resu.setname(name)
                results.append(resu)
            multiple_cor(results,find_output(key,'bc2', out_path))    
            plot_histoMU(results,nbreads,find_output(key,'bc1', out_path))
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
            nbreads=calcnbread(data_path[len(data_path)-1],key)
            for tool in lst:
                save = pysam.set_verbosity(0)
                bamFP = pysam.AlignmentFile(tool, "rb")
                pysam.set_verbosity(save)
                name=findname(tool)
                resu=countdiff(bamFP)
                resu.setname(name)
                results.append(resu)
            common_error_gp1(results,find_output(key,'rl4', out_path))
            common_error_gp2(results,find_output(key,'rl5', out_path))
            multiple_cor(results,find_output(key,'rl2', out_path))    
            plot_histoMU(results,nbreads,find_output(key,'rl1', out_path))
            plot_histoNotHuman(results,find_output(key,'rl3', out_path))
            if test:
                return results

def calcnbread(files,key):
    pos=list(open(files).readlines())
    for line in pos:
        if key in line:
            return int(line[line.index(':')+1:])            

