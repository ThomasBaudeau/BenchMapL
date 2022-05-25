from distutils.command.config import config
import pysam
import numpy as np
import matplotlib.pyplot as plt
from BenchPlot import *

def parsepath(path):
    one=(path.index("/"))
    two=(path.index('_',one))
    return path[one+1:two]           

def group_input(files):
    group={}
    for specie in snakemake.config["species"]:
        for length in snakemake.config['length']:
            for er in snakemake.config['error_rate']:
                for file in files:
                    if  file.count('_'+specie+'_')==1 and file.count('_'+str(length)+'_')==1  and file.count('_'+str(er)+'_')==1:
                        try:
                            group[specie+'_'+str(length)+'_'+str(er)+'_'].append(file)
                        except:
                            group[specie+'_'+str(length)+'_'+str(er)+'_']=[]
                            group[specie+'_'+str(length)+'_'+str(er)+'_'].append(file)
    return group

def do_something(data_path, out_path, myparam):
    files=list(snakemake.input)
    groups=group_input(files)
    for key in groups.keys():
        if find_output(key,outpath=out_path): 
            lst=groups[key]
            results=[]
            for tool in lst:
                save = pysam.set_verbosity(0)
                bamFP = pysam.AlignmentFile(tool, "rb")
                pysam.set_verbosity(save)
                name=parsepath(tool)
                resu=countdiff(bamFP)
                resu.setname(name)
                results.append(resu)
            common_error_gp1(results,find_output(key,'rl4', out_path))
            common_error_gp2(results,find_output(key,'rl5', out_path))
            multiple_cor(results,find_output(key,'rl2', out_path))    
            plot_histoMU(results,find_output(key,'rl1', out_path))
            plot_histoNotHuman(results,find_output(key,'rl3', out_path))

do_something(snakemake.input, snakemake.output, snakemake.config["param"])