import os
import shutil
import glob

os.chdir('Workflow')
try:
    shutil.rmtree('result')
except:
    pass
print("1%")
try:
    shutil.rmtree('perfect_sam')
except:
    pass
try:
    shutil.rmtree('medaka')
except:
    pass
try:
    shutil.rmtree('perfect_sam')
except:
    pass

try:
    shutil.rmtree('mapped_reads')
except:
    pass
print("80%")
try:
   os.remove('variant_file.txt')
except:
    pass

try:
   os.remove('data/variant_file.txt')
except:
    pass

try:
   os.remove()
except:
    pass
try:
    shutil.rmtree('data/samples')
except:
    pass
try:
    shutil.rmtree('data/perfect')
except:
    pass
for f in glob.glob('data/*.fai'):
    os.remove(f)
for f in glob.glob('data/ref_*_*'):
    os.remove(f)

for f in glob.glob('*_k10.txt'):
    os.remove(f)

try:
   os.remove('result.csv')
except:
    pass
