from distutils.command.config import config
import pysam
import numpy as np
import matplotlib.pyplot as plt
from BenchPlot import *

def parsepath(path):
    """check the command and the tools

    :param path: a path
    :type path: string
    :return: the string name of the command
    :rtype: string
    """
    asplit=path.split('_')
    command=asplit[5][:-4]
    tool=asplit[1][6:]
    for pos in range(len(snakemake.config['param'][tool]['command'])):
        if snakemake.config['param'][tool]['command'][pos]==command:
            return snakemake.config['param'][tool]['name'][pos]
    raise ValueError("Wrong Input")
           

def group_input(files):
    """generate a dict for each group of tools/datasets/command

    :param files: list of bam file from the different tools/dataset
    :type files: list[string]
    :return: dictionnary of the different datasets groups with each files in a list
    :rtype: dict
    """
    group={}
    for specie in snakemake.config["species"]:
        for length in snakemake.config['length']:
             for er in snakemake.config['error_rate']:
                for t in snakemake.config['param']:
                    for file in files:
                        if  file.count('_'+specie+'_')==1 and file.count('_'+str(length)+'_')==1  and file.count('_'+str(er)+'_')==1 and file.count(t+'_')==1:
                            try:
                                group[t+'_'+specie+'_'+str(length)+'_'+str(er)+'_'].append(file)
                            except:
                                group[t+'_'+specie+'_'+str(length)+'_'+str(er)+'_']=[]
                                group[t+'_'+specie+'_'+str(length)+'_'+str(er)+'_'].append(file)
    return group

def do_something(data_path, out_path, myparam):
    """main function

    :param data_path: snakemake.input
    :type data_path: list
    :param out_path: snakemake.output
    :type out_path: list
    :param myparam: config files
    :type myparam: list
    """
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
            multiple_cor(results,find_output(key,'bc2', out_path))    
            plot_histoMU(results,find_output(key,'bc1', out_path))
            common_error_gp2(results,find_output(key,'bc3', out_path))

do_something(snakemake.input, snakemake.output, snakemake.config["param"])
