#BSUB-J MOLM13_alignment               # name of the job
#BSUB-W 3:00                   	# Time limit in minutes
#BSUB-n 4
#BSUB-R "span[ptile=4]"
#BSUB-M 10

#BSUB-o stdout
#BSUB-o stderr

# STAR alignment
ls -d Sample*[0-9] | xargs -I {} sh -c "STAR --genomeLoad NoSharedMemory --genomeDir hg19_star/ --readFilesIn {}/*R1_001.fastq.gz {}/*R2_001.fastq.gz --runThreadN 4 --alignIntronMin 70 --alignIntronMax 100000 --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 1 --outFilterMultimapScoreRange 0 --outFilterMismatchNmax 5 --outFileNamePrefix {}/ --outStd Log --readFilesCommand zcat"
ls -d Sample*[0-9] | xargs -I {} sh -c "mv {}/*bam {}.bam"
ls Sample*[0-9].bam | xargs -n 1 samtools index
