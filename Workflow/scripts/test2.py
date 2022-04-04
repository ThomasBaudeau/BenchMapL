from distutils.command.config import config
import pysam
import numpy as np
import matplotlib.pyplot as plt
Tools=['mapped_reads/graphmap_default_align#-x#sensitive.bam','mapped_reads/mm2f_default_-ax#map-ont.bam','mapped_reads/minimap2_default_-ax#map-ont.bam','mapped_reads/graphmap2_default_align#-a#sggotoh.bam','mapped_reads/blasr_default_--sam.bam']

class result:
    def __init__(self):
       
        self.cor=0
        self.cor_5=0
        self.cor_10=0
        self.cor_20=0
        self.mapped=0
        self.unmapped=0

    def setname(self,name):
         self.name=name

    def corpercent(self,cor):
        return (cor/(self.mapped+self.unmapped))*100

    def mappercentU(self):
        return (self.unmapped/(self.mapped+self.unmapped))*100

    def mappercentM(self):
        return (self.mapped/(self.mapped+self.unmapped))*100


    def increm_cor(self,num=0):
        if num not in [0,5,10,20]:
            print("Incorect Input")
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



def parsepath(path):
    one=(path.index("/"))
    two=(path.index('_',one))
    return path[one+1:two]
    

def parsenamesimu(text):
    name=text
    st_pos= int(text.split('_')[1])#int(text[text.index('_')+1:text.index(';')])
    end_pos=int(text.split('_')[6])#int(text[text.index(';')+1:len(text)-1].split('_')[4])

    return st_pos,end_pos,(st_pos+end_pos)

def read_align(resu,read,threshold=0):
    st_ref,mbase_ref,end_ref=parsenamesimu(read.query_name)
    st_alg=read.reference_start
    end_alg=read.reference_end
    mbase_alg=read.reference_length
    crit1=abs(st_ref-st_alg)
    crit2=abs(end_ref-end_alg)
    crit3=abs(mbase_ref-mbase_alg)
    if crit1<=threshold and crit2<=threshold and crit3<=threshold:
        resu.increm_cor(threshold)
        return True
    else:
        return False

def countdiff(bamFP):
    resu=result()
    for read in bamFP:
        if not read.is_unmapped:
            resu.increm_mapped()
            for threshold in [0,5,10,20]:
                rep=read_align(resu,read,threshold)
                if rep:
                    break
        else:
            resu.increm_unmapped()
    return resu

def find_output(key,name):
    print(snakemake.output)
    print(key)
    for output in snakemake.output:
        if output.count(key)!=0 and output.count(name)!=0:
            return output
    return False
            

def plot_histoMU(results,key):
    output=find_output(key,'test1')
    if output:
        width = 0.35       # the width of the bars: can also be len(x) sequence
        fig, ax = plt.subplots()
        labels=get_label(results)
        d1,d2=get_MapOrNot(results)
        ax.bar(labels, d2, width, label='Mapped_reads')
        ax.bar(labels, d1, width, label='Unmapped_reads',bottom=d2)
        ax.set_ylabel('nb_reads')
        ax.set_title('number of reads by status and tools')
        ax.legend()
        plt.savefig(output)

def get_label(results):
    labels=[]
    for result in results:
        labels.append(result.name)
    return labels

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


def multiple_cor(results,key):
    output=find_output(key,'test2')
    if output:
        data = get_all_cor(results)

        columns = get_label(results)
        rows = ['perfect','5pb_errors','10pb_errors','20pb_errors']

        values = np.arange(0, 100, 5)

        # Get some pastel shades for the colors
        colors = plt.cm.BuPu(np.linspace(0.15, .65, len(rows)))
        n_rows = len(data)

        index = np.arange(len(columns)) + 0.3
        bar_width = 0.4

        # Initialize the vertical-offset for the stacked bar chart.
        y_offset = np.zeros(len(columns))
        y=np.zeros(len(columns))
        # Plot bars and create text labels for the table
        cell_text = []
        for row in range(n_rows):
            plt.bar(index, data[row], bar_width, bottom=y, color=colors[row])
            y_offset = data[row]
            y=y+data[row]
            cell_text.append(['%1i' % (x / 1) for x in y])
        # Reverse colors and text labels to display the last value at the top.



        # Add a table at the bottom of the axes
        the_table = plt.table(cellText=cell_text,
                            rowLabels=rows,
                            rowColours=colors,
                            colLabels=columns,
                            loc='bottom')

        # Adjust layout to make room for the table:
        plt.subplots_adjust(left=0.2, bottom=0.2)

        plt.ylabel("nb reads correctly aligned")
        plt.yticks(values, ['%d' % val for val in values])
        plt.xticks([])
        plt.title('proportion of correctly mapped read by tools and error')

        plt.savefig(output)

def group_input(files):
    group={}
    for specie in snakemake.config["species"]:
        for length in snakemake.config['length']:
            for file in files:
                if  file.count('_'+specie+'_')==1 and file.count('_'+str(length)+'_')==1 :
                    try:
                        group[specie+'_'+str(length)+'_'].append(file)
                    except:
                        group[specie+'_'+str(length)+'_']=[]
                        group[specie+'_'+str(length)+'_'].append(file)
    return group

def do_something(data_path, out_path, myparam):
    files=list(snakemake.input)
    groups=group_input(files)
    for key in groups.keys(): 
        lst=groups[key]
        results=[]
        for tool in lst:
            bamFP = pysam.AlignmentFile(tool, "rb")
            name=parsepath(tool)
            count=0
            resu=countdiff(bamFP)
            resu.setname(name)
            results.append(resu)
        multiple_cor(results,key)    
        plot_histoMU(results,key)

do_something(snakemake.input[0], snakemake.output[0], snakemake.config["param"])


