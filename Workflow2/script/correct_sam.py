def find_variant():
    vf=open('variant_file.txt').readlines()
    vardict={}
    for line in vf:
        lsplit=line.split('\t')
        if lsplit[0] not in vardict.keys():
            vardict[lsplit[0]]=[[],[]]
    pass