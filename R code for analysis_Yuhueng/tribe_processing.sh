#BSUB-J rnaseq_alignment               # name of the job
#BSUB-W 2:30                   	# Time limit in minutes
#BSUB-n 6
#BSUB-R "span[ptile=6]"
#BSUB-M 6

#BSUB-o stdout
#BSUB-o stderr

# STAR alignment
ls -d Sample* | xargs -I {} sh -c "STAR --genomeLoad NoSharedMemory --genomeDir mm10_star/ --readFilesIn {}/*R1_001.fastq.gz {}/*R2_001.fastq.gz --runThreadN 6 --alignIntronMin 70 --alignIntronMax 100000 --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 1 --outFilterMultimapScoreRange 0 --outFilterMismatchNmax 5 --outFileNamePrefix {}/ --outStd Log --readFilesCommand zcat"
ls -d Sample* | xargs -I {} sh -c "mv {}/*bam {}.bam"
ls Sample*bam | xargs -n 1 samtools index
