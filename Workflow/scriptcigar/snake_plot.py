from plotcigar import *

def do_something(input,output,params):
    print(output)
    mean_score(input[0:len(input)-2],input[-1],output,base=params)


do_something(snakemake.input, snakemake.output, snakemake.params)