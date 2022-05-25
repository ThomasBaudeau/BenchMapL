from distutils.command.config import config
import pysam
import numpy as np
import matplotlib.pyplot as plt
from BenchPlot import *
       

def group_input(files):
    group={}
    for specie in ['vih']:
        for length in ['350']:
            for er in ['0.9']:
                for file in files:
                    if  file.count('_'+specie+'_')==1 and file.count('_'+str(length)+'_')==1  and file.count('_'+str(er)+'_')==1:
                        try:
                            group[specie+'_'+str(length)+'_'+str(er)+'_'].append(file)
                        except:
                            group[specie+'_'+str(length)+'_'+str(er)+'_']=[]
                            group[specie+'_'+str(length)+'_'+str(er)+'_'].append(file)
    return group

def do_something(data_path, out_path, myparam):
    files=[data_path]
    print(files)
    groups=group_input(files)
    print(groups)
    for key in groups.keys():
            lst=groups[key]
            print(lst)
            results=[]
            for tool in lst:
                save = pysam.set_verbosity(0)
                bamFP = pysam.AlignmentFile(tool, "rb",check_sq=False)
                pysam.set_verbosity(save)
                name='outil'
                resu=countdiff(bamFP)
                resu.setname(name)
                results.append(resu)
            print(results)
            common_error_gp1(results,find_output(key,'rl4', out_path))
            common_error_gp2(results,find_output(key,'rl5', out_path))
            multiple_cor(results,find_output(key,'rl2', out_path))    
            plot_histoMU(results,find_output(key,'rl1', out_path))
            plot_histoNotHuman(results,find_output(key,'rl3', out_path))

do_something('tool_vih_350_0.9_FF.sam', 'okresu.png', 'ok')