#!/bin/bash

input_file=$(realpath $1)
output_file=$(realpath $2)

awk '
BEGIN { FS=OFS="\t" } 
($1 == "#CHROM") {a=1}
a==1 {$10=$NF; for (i = 11; i <= NF; i++) {$i=""}}
{print}' "$input_file" | sed -r 's/QNAME=[a-z A-Z 0-9]*[@][0-9]+\t*/QNAME=species@1/'|sed -r 's/[0-9]+:sample\t*/1:sample/'|sed -r 's/.\/.\t*/1/'> "$output_file"
