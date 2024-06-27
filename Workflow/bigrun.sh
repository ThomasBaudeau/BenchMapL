#!/bin/bash

var1=../tool/BenchMapL/config
var2=../../../result
count=0
tt= ls |  wc -l 

conda activate snakemake;
for f in *;
	do cp $f $var1/;
	count=$(($count+1))
	mv $var1/$f $var1/config.yaml ;
	cd ../tool/BenchMapL/Workflow;
	if ! snakemake -c8 all --use-conda; then
        echo "Faut tout relancer" | mail -r bws@univ-lille.fr thomas.baudeau@univ-lille.fr -s "Plantage" ;
        exit 0;
    fi
	mkdir $var2/config$count;
    cp $var1/config.yaml $var2/config$count/config.yaml
	mv resultfinal $var2/config$count/;
	mv bcf $var2/config$count/
    mv medaka $var2/config$count/
	mv plots $var2/config$count/
	mv benchmarks $var2/config$count/
	mv mapped_reads $var2/config$count/
	mv data/samples/ $var2/config$count/
    python3 clear.py;
	cd ../../../data;
	echo echo $count / $tt
done ||exit 0
echo "FIN" | mail -r bws@univ-lille.fr thomas.baudeau@univ-lille.fr -s "FIN" ;