import re
import matplotlib.pyplot as plt
import pysam
import numpy as np
from upsetplot import from_contents
from upsetplot import UpSet

class result:
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
         self.name=name

    def corpercent(self,cor):
        return 0 if self.mapped==0 else (cor/(self.mapped))*100

    def mappercentU(self):
        return (self.unmapped/(self.mapped+self.unmapped))*100

    def mappercentM(self):
        return (self.mapped/(self.mapped+self.unmapped))*100

    def addGroup1(self,add):
        self.group1.append(add)

    def addGroup2(self,add):
        self.group2.append(add)

    def missalign(self):
        self.missaligned+=1

    def increm_cor(self,num=0):
        if num not in [0,5,10,20]:
            print("Wrong Input")
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
        self.mapped+=1

    def increm_unmapped(self):
        self.unmapped+=1

def parsenamesimu(text):
    st_pos= int(text.split('_')[1])
    end_pos=int(text.split('_')[6])
    return st_pos,end_pos,(st_pos+end_pos)

def read_align(resu,read,threshold=0):
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
    for output in outpath:
        if name:
            if output.count(key)!=0 and output.count(name)!=0:
                return output
        else:
            if output.count(key)!=0:
                return True
    return False

def get_label(results):
    labels=[]
    for result in results:
        labels.append(result.name)
    return labels

def is_human(read,resu):
    if read.query_name.split('_')[2]=='human':
        if not read.is_unmapped:
            resu.missalign()
        return False
    return True

def get_nbHuman(results):
    nb=[]
    for result in results:  
        nb.append(result.missaligned)
    return nb

def get_MapOrNot(results):
    data_unmap=[]
    data_map=[]
    for result in results:
        data_unmap.append(result.mappercentU())
        data_map.append(result.mappercentM())
    return data_unmap,data_map

def get_all_cor(results):
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

def plot_histoMU(results,output):
    if output:
        width = 0.45       # the width of the bars: can also be len(x) sequence
        fig, ax = plt.subplots(figsize=(7, 7))
        plt.rcParams.update({'font.size': 11})
        labels=get_label(results)
        d1,d2=get_MapOrNot(results)
        ax.set_xticklabels(labels, rotation = 45)
        ax.bar(labels, d2, width, label='Mapped_reads')
        ax.bar(labels, d1, width, label='Unmapped_reads',bottom=d2)
        ax.set_ylabel('Percentage of reads')
        ax.set_title('Percentage of reads by status and tools')
        ax.legend()
        plt.subplots_adjust(left=0.2, bottom=0.2)
        plt.savefig(output,dpi=300,format='pdf')

def plot_histoNotHuman(results,output):
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