import os
import shutil
import glob

def do_something(data_path, out_path,param):
    os.renames(data_path[0],out_path[0])
    os.renames(data_path[0].replace('variantreads.fastq','perfect.sam'),out_path[1])
    shutil.rmtree('result')
    shutil.rmtree('perfect_sam')


do_something(snakemake.input,snakemake.output,snakemake.config)