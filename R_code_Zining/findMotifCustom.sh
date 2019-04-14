#!/bin/bash

# input information
inputpar=./fasta
inputfilepath=mpp2_cds_target.fa
bgfilepath=mpp2_cds_background.fa
input="$inputpar/$inputfilepath"
bg="$inputpar/$bgfilepath"

# output information
par=./mpp2_cds/
denovodirectory=denovo
knowndirectory=known

# motif information
motif=./homer/motifs/cisbp.motif

./homer/bin/findMotifs.pl $input fasta $par$denovodirectory -fasta $bg -noknown -len 6,7,8 -rna -norevopp -nlen 2 -olen 2 
./homer/bin/compareMotifs.pl $par$denovodirectory$"/homerMotifs.all.motifs" $par -rna -known $motif -cpu 4 -norevopp
# ./homer/bin/findMotifs.pl $motiffasta fasta $par -norevopp -find $par$denovodirectory$"/homerResults/motif1.motif" > $par$"all_found.txt"
