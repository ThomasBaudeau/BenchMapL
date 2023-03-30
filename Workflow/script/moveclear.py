import os
import shutil
import glob

def do_something(data_path, out_path,param):
    todel='data/'+data_path[0].split('/')[2].split('_')[0]+'_'
    os.renames(data_path[0],out_path[0])
    os.renames(data_path[0].replace('variantreads.fastq','perfect.sam'),out_path[1])
    shutil.rmtree('result')
    shutil.rmtree('perfect_sam')
    os.renames('variant_file.txt','data/variant_file.txt')
    for f in glob.glob(todel+'*'):
        os.remove(f)


do_something(snakemake.input,snakemake.output,snakemake.config)