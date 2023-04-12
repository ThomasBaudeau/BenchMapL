from plotcigar import *

def do_something(input,output,params):
    onexp={}
    for abam in input:
        firstpart=abam.split("mapped_reads/")[1]
        secondpart=firstpart.split('_')
        xp='_'.join(secondpart[1:4])
        if secondpart[0]=='perfect':
            try:
                onexp[xp]['perfect']=abam
            except:
                onexp[xp]={}
                onexp[xp]['perfect']=abam
        else:
            try:
                onexp[xp]['lstbam'].append(abam)
            except:
                onexp[xp]={}
                onexp[xp]['lstbam']=[abam]
    
    overall(onexp,str(output))
 

do_something(snakemake.input, snakemake.output, snakemake.params)
