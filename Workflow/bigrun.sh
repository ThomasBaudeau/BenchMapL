#!/bin/bash

var1=../TODO
var2=../TODO
count=0
tt= ls |  wc -l 

conda activate snakemake;
for f in *;
	do cp $f $var1;
	count=$(($count+1))
	mv $var1/$f $var1/config.yaml ;
	cd TODO;
	snakemake -c8 all --use-conda ||(echo "Faut tout relancer" | mail -r bws@univ-lille.fr thomas.baudeau@univ-lille.fr -s "Plantage" ; exit 0);
    cd TODO
    mkdir /TODO/config$count;
    cp $var1/config.yaml TODO/config$count/config.yaml
	mv resultfinal TODO/config$count/;
	mv bcf TODO/config$count/
    mv medaka TODO/config$count/
    python3 clear.py;
	cd $var2;
	echo echo $count / $tt
done ||exit 0
echo "FIN" | mail -r bws@univ-lille.fr thomas.baudeau@univ-lille.fr -s "FIN" ;