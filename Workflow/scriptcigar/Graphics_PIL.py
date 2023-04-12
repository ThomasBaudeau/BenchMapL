from PIL import Image, ImageDraw, ImageFont
from extract import scores
from math import *
import pandas as pd
import matplotlib.pyplot as plt

# # Files to use
# file1 = "Scripts/perfect_vih_350_0.9_perfect.bam"
# file2 = "Scripts/minimap2_vih_350_0.9_-k8-w1.bam"
# file3 = "Scripts/graphmap2_vih_350_0.9_.bam"
# file4 = "Scripts/magicblast_vih_350_0.9_-spliceF-gapopen4-gapextend2.bam"


# Extraction and score implementation
# l_nucl = scores(file2, perfect=file1)
#Coloriage
def color_car(score):
    if score == 1:
       return '#6aff6a'
    elif score == 0:
       return '#ff6961'
    else:
       return '#ffffff'

# Contours et délimitation des carrés de chaque colonne


def find_max_cov(lnucl):
    max_rec=0
    for i in lnucl:
        rec=len(i.get_score())
        if rec > max_rec:
            max_rec=rec
    return max_rec+5



def visualisation(l_nucl, outfile,height = 15000 ):
    cov=find_max_cov(l_nucl)
    decoupe= int(sqrt(((len(l_nucl)+5)*5)/((cov+5)*5)))+1
    height=int((decoupe*(cov+5)*5)+50)+1
    width = int((len(l_nucl)*5)/decoupe)+150
    size = (width, height)
    image = Image.new('RGB', size,'white')
    draw = ImageDraw.Draw(image)
    font = ImageFont.truetype('Scripts/verdana.ttf', 8)
    c = 50
    h=1
    px = c
    nb_rec = 1
    lh=0
    max_rec=0
    d=c+lh
    for i in l_nucl:
        if i.get_pos()%5 ==0 :
            draw.line([(px+3,height-(d-2)), (px+3,height-(d-8))], fill='black')
            draw.text((px, height-(d-10)), f"{i.get_pos()}", font=font, fill ="black", align="center")
        
        for rec in i.get_score():
            draw.rectangle([px, (height-d)-nb_rec*5, px+5, (height-d)-(nb_rec-1)*5], fill=color_car(rec), width=1)
            nb_rec += 1
        draw.rectangle([px, (height-d)-i.get_count()*5, px+5, (height-d)], fill=None, outline='black', width=1)
        if nb_rec>max_rec:
            max_rec=nb_rec
        nb_rec = 1
        px+=5
        if (i.get_pos())%(int((width-50)/5)-10) == 0:
            px = c
            lh=5*(max_rec+5)
            d=d+c+lh
            # image.save('Comparison_Score&Coverage_Tools_to_pos{}.png'.format(i.get_pos()))
            # image.show()
            # image = Image.new('RGB', size,'white')
            # draw = ImageDraw.Draw(image)
            
        elif i.get_pos() == (len(l_nucl)-1):
            image.save(outfile)

    

def meanpos(lnucl):
    mscore=[]
    for i in lnucl:
        mscore.append(i.get_mean_score())
    return mscore


def meanposlisse(lnucl,lissage=10):
    mscore=[]
    tt=0
    ttlen=0
    for idx,i in enumerate(lnucl):
        tt+=sum(i.get_score())
        ttlen+=len(i.get_score())
        if idx%lissage==0:
            if ttlen>0:
                amean=tt/ttlen

            else:
                amean=0    
            mscore.append(amean)
            tt=0
            ttlen=0
    return mscore


def meanposflot(lnucl,lissage=10):
    mscore=[]
    flotlist=[]
    tt=0
    ttlen=0
    for idx,i in enumerate(lnucl):
        lasomme=sum(i.get_score())
        lalong=len(i.get_score())
        tt+=lasomme
        ttlen+=lalong
        flotlist.append((lasomme,lalong))
        if idx>9:
            a=flotlist.pop(0)
            tt-=a[0]
            ttlen-=a[1]
        if ttlen>0:
            amean=tt/ttlen
        else:
            amean=0    
        mscore.append(amean)
    return mscore

def graph_score(dic,outfile):
    df=pd.DataFrame.from_dict(dic)
    df.plot(subplots=True, linewidth=0.5)
    plt.savefig(str(outfile))

def overall_score(lnucl):
    ttsum=0
    ttlen=0
    for nucl in lnucl:
        ttsum+=sum(nucl.get_score())
        ttlen+=len(nucl.get_score())
    if ttlen>0:
        return ttsum/ttlen
    else:
        return 0
    





