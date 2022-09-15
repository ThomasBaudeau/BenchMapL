import re
import matplotlib.pyplot as plt
import pysam
from upsetplot import from_contents
from upsetplot import UpSet
import os
import sys
sys.path.insert(1, '../Workflow/scripts')
import BenchPlot


CONFIG_FILE={'length':[5000],
'species':['test'],
'error_rate':['5'],
'param:':[{'minimap':{'name':['default'],'command':['#']}},
{'minimap2':{'name':['default'],'command':['#']}},
{'graph':{'name':['default'],'command':['#']}},
{'graph2':{'name':['default'],'command':['#']}},
{'blost':{'name':['default'],'command':['#']}}]}

EXPECTED_RESULT= {'minimap':[100,0,0,0,100,0,0,0,0],'minimap2':[90,0,0,0,90,10,0,0,0],'graph':[80,0,0,0,80,20,0,0,0],'graph2':[70,0,0,0,70,30,0,0,0],'blost':[60,0,0,0,60,40,0,0,0]}
EXPECTED_RESULT2= {'minimap':[70,10,10,10,100,0,0,0,0],'minimap2':[30,20,20,10,90,10,10,0,0],'graph':[20,20,20,20,80,20,0,0,0],'graph2':[25,10,15,20,70,30,0,0,0],'blost':[60,0,0,0,60,40,0,0,0]}
EXPECTED_RESULT3= {'minimap':[70,10,10,10,100,0,0,0,0],'minimap2':[30,20,20,10,90,10,10,10,0],'graph':[20,20,20,20,80,20,0,0,20],'graph2':[25,10,15,20,70,30,0,10,0],'blost':[60,0,0,0,60,40,0,10,10]}



def build_sam(a,name,flag,start,cigar=((0,5),(0,5))):
        a.query_name = "read"+name
        a.query_sequence="AAAAAAAAAA"
        a.flag = flag
        a.reference_id = 0
        a.reference_start = start
        a.mapping_quality = 2
        a.cigar = cigar
        a.next_reference_id = 0
        a.next_reference_start=0
        a.template_length=167
        a.query_qualities = pysam.qualitystring_to_array("<<<<<<<<<<")
        a.tags = (("NM", 1),
                ("RG", "L1"))
        return a
def build_files(size,command,species,er,expected):
    us='_'
    lst=[]
    header = { 'HD': {'VN': '1.0'},
            'SQ': [{'LN': 1575, 'SN': 'chr1'},
                   {'LN': 1584, 'SN': 'chr2'}] }
    for tool in EXPECTED_RESULT.keys():
        filename=tool+us+species+us+size+us+er+us+command+'.sam'
        with pysam.AlignmentFile('READTEST/'+filename, "w", header=header) as outf:
            a = pysam.AlignedSegment()
            for i in range(expected[tool][0]):
                a=build_sam(a,'a'+str(i)+'_'+'10_aligned_-1_F_-1_10_20_c33',0,10)
                outf.write(a)
            for j in range(expected[tool][1]):
                a=build_sam(a,'b'+str(j)+'_'+'10_aligned_-1_F_-1_10_20_c33',0,12)
                outf.write(a)
            for k in range(expected[tool][2]):
                a=build_sam(a,'c'+str(k)+'_'+'10_aligned_-1_F_-1_10_20_c33',0,14)
                outf.write(a)
            for l in range(expected[tool][3]):
                a=build_sam(a,'d'+str(l)+'_'+'10_aligned_-1_F_-1_10_20_c33',0,18)
                outf.write(a)
            for m in range(expected[tool][5]):
                a=build_sam(a,'e'+str(m)+'_'+'-1000_aligned_-1_F_10_-1_-2000_c33',4,1000)
                outf.write(a)
            for n in range(expected[tool][6]):
                a=build_sam(a,'f'+str(n)+'_'+'10_aligned_-1_F_-1_10_20_c33',0,50)
                outf.write(a)
            for o in range(expected[tool][8]):
                a=build_sam(a,'e'+str(o)+'_'+'-1000_human_-1_F_10_-1_-2000_c33',4,1000)
                outf.write(a)
            for p in range(expected[tool][7]):
                a=build_sam(a,'b'+str(p)+'_'+'10_human_-1_F_-1_10_20_c33',0,12)
                outf.write(a)
        lst.append('READTEST/'+filename)
    return lst


def test1():
    outpath=['test_5000_5_rl1.pdf','test_5000_5_rl2.pdf','test_5000_5_rl3.pdf','test_5000_5_rl4.pdf','test_5000_5_rl5.pdf']
    apath=build_files('5000','#','test','5',EXPECTED_RESULT)
    ok=BenchPlot.plot_simple(apath,outpath,CONFIG_FILE,True)
    for obj in ok:
        assert EXPECTED_RESULT[obj.name][0]==obj.mapped
        assert EXPECTED_RESULT[obj.name][5]==obj.unmapped
    for file in apath:
        os.remove(file)
    for file in outpath:
        os.remove(file)


def test2():
    outpath=['test_5000_5_rl1.pdf','test_5000_5_rl2.pdf','test_5000_5_rl3.pdf','test_5000_5_rl4.pdf','test_5000_5_rl5.pdf']
    apath=build_files('5000','#','test','5',EXPECTED_RESULT2)
    ok=BenchPlot.plot_simple(apath,outpath,CONFIG_FILE,True)
    for obj in ok:
        assert EXPECTED_RESULT2[obj.name][0]==obj.cor
        assert EXPECTED_RESULT2[obj.name][1]==obj.cor_5,(obj.cor,obj.cor_5,obj.cor_10,obj.cor_20)
        assert EXPECTED_RESULT2[obj.name][2]==obj.cor_10,(obj.cor,obj.cor_5,obj.cor_10,obj.cor_20)
        assert EXPECTED_RESULT2[obj.name][3]==obj.cor_20,(obj.cor,obj.cor_5,obj.cor_10,obj.cor_20)
        assert EXPECTED_RESULT2[obj.name][6]==(obj.mapped-(obj.cor+obj.cor_5+obj.cor_10+obj.cor_20))
        assert EXPECTED_RESULT2[obj.name][4]==obj.mapped,(obj.mapped,obj.name,EXPECTED_RESULT2[obj.name][4])
        assert EXPECTED_RESULT2[obj.name][5]==obj.unmapped
    for file in apath:
        os.remove(file)
    for file in outpath:
        os.remove(file)

def test3():
    outpath=['test_5000_5_rl1.pdf','test_5000_5_rl2.pdf','test_5000_5_rl3.pdf','test_5000_5_rl4.pdf','test_5000_5_rl5.pdf']
    apath=build_files('5000','#','test','5',EXPECTED_RESULT3)
    ok=BenchPlot.plot_simple(apath,outpath,CONFIG_FILE,True)
    for obj in ok:
        assert EXPECTED_RESULT3[obj.name][0]==obj.cor
        assert EXPECTED_RESULT3[obj.name][1]==obj.cor_5,(obj.cor,obj.cor_5,obj.cor_10,obj.cor_20)
        assert EXPECTED_RESULT3[obj.name][2]==obj.cor_10,(obj.cor,obj.cor_5,obj.cor_10,obj.cor_20)
        assert EXPECTED_RESULT3[obj.name][3]==obj.cor_20,(obj.cor,obj.cor_5,obj.cor_10,obj.cor_20)
        assert EXPECTED_RESULT3[obj.name][6]==(obj.mapped-(obj.cor+obj.cor_5+obj.cor_10+obj.cor_20))
        assert EXPECTED_RESULT3[obj.name][4]==obj.mapped,(obj.mapped,obj.name,EXPECTED_RESULT3[obj.name][4])
        assert EXPECTED_RESULT3[obj.name][5]==obj.unmapped
        assert EXPECTED_RESULT3[obj.name][7]==obj.missaligned
    for file in apath:
        os.remove(file)
    for file in outpath:
        pass

test1()
test2()
test3()