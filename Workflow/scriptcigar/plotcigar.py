from extract import scores
from tqdm import tqdm
from Graphics_PIL import visualisation
from Graphics_PIL import meanpos
from Graphics_PIL import meanposlisse
from Graphics_PIL import meanposflot
from Graphics_PIL import graph_score
from Graphics_PIL import overall_score
import numpy as np


def plot_histcig(f,perfect,outfile):
    # parser=argparse.ArgumentParser(prog = "__main__.py", description = "Visualisation of a comparaison between a theorical perfect alignment and a alignement coming from a mapping tools")

    # # Required arguments
    # parser.add_argument("BAMfiles",nargs='*', help = "path to one or several file in BAM format (.bam) containing the aligend reads")
    # parser.add_argument("-perfect", help = "path to the  file in BAM format (.bam) containing the theorical perfect alignement i", required = True)
    # parser.add_argument("-heigth", help = "heigth of the outcome images", type = int, default=300)
    # args = parser.parse_args(sys.argv[1:])

    # Display
    l_nucl = scores((f,), perfect)
       

    visualisation(l_nucl,outfile)

# import doctest
# doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE | doctest.ELLIPSIS, verbose=True)
    
# if __name__ == "__main__":
#     main()




def mean_score(ltools,perfect,outfile,base='#'):
    sc={}
    ov={}
    for i in ltools:
        l_nucl = scores((i,), perfect)
        ov[i.split('mapped_reads/')[1]]=overall_score(l_nucl)
        base=str(base)
        a=i.split('mapped_reads/')[1].split('_')[4].replace('.bam', '')
        if a == base:
            nor = meanposflot(l_nucl)
        else:
            sc[i.split('mapped_reads/')[1]]=meanposflot(l_nucl)
    for i in sc.keys():
        sc[i]=(np.array(sc[i])-np.array(nor)).tolist()    
    graph_score(sc,outfile)


def overall(adic,outfile): 
    file=open(outfile,'w')
    ov={}
    for xp in tqdm(adic.keys()):
        perfect= adic[xp]['perfect']
        for i in adic[xp]['lstbam']:
            l_nucl = scores((i,), perfect)
            ov[i.split('mapped_reads/')[1]]=overall_score(l_nucl)
        for x in ov.keys():
            file.write(str(x)+"\t"+str(ov[x])+'\n')
