#!/bin/bash

## install fastp and use fastp filter data and merge R1 R2 data
dir=
name=$1

date

fastp  -i ${dir}/${name}/*${name}_combined_R1.fastq.gz   -I ${dir}/${name}/*${name}_combined_R2.fastq.gz   --merge --merged_out ${dir}/${name}/${name}_merged_out.fastq.gz --include_unmerged ${dir}/${name}/${name}_unmerged_out.fastq.gz  -j ${dir}/${name}/${name}.json -h ${dir}/${name}/${name}.html -z 9 -f 0 -t 0 -T 0 -F 0 -L   --thread=10

date
echo "Done"
