import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def parseSimplestats(file):
    data=[]
    files= open(file)
    lignes=files.readlines()
    lignes=lignes[0:50]
    lignes=list(filter(check_SN,lignes))
    lignes=[l.split('\t' )for l in lignes]
    print(lignes)
    for l in lignes:
       a=filter(is_numeric,l)
       for i in a:
           data.append(int(i.replace('\n','')))
    print(data)
    return data

def check_SN(n):
    if 'SN'==n[0:2]:
        return n
def is_numeric(n):
    num=n.replace('\n','')
    if num.isnumeric() or num=='0':
       return True 


def parsepath(path):
    print(path)
    one=(path.index("/"))
    two=(path.index('_',one))
    return path[one+1:two]
  




def plot_histo(labels,data):
    print(data)
    width = 0.35       # the width of the bars: can also be len(x) sequence
    fig, ax = plt.subplots()
    print(data)
    ax.bar(labels, data, width, label='Unmapped_reads')
    ax.set_ylabel('nb_reads')
    ax.set_title('Scores by group and gender')
    ax.legend()

    plt.show()
    plt.save()


def do_something(data_path, out_path, myparam):
    files=list(snakemake.input)
    data=[]
    label=[]
    for path in files:
        data.append(parseSimplestats(path)[8])
        label.append(parsepath(path))
    plot_histo(label,data)

do_something(snakemake.input[0], snakemake.output[0], snakemake.config["param"])