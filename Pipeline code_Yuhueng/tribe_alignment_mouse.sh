#BSUB-J LSC_alignment               # name of the job
#BSUB-W 2:00                   	# Time limit in minutes
#BSUB-n 9
#BSUB-M 4
#BSUB-R "span[ptile=9]"

#BSUB-o stdout
#BSUB-o stderr

# STAR alignment
ls -d Sample_*[0-9] | xargs -I {} sh -c "STAR --genomeLoad NoSharedMemory --genomeDir mm10_star/ --readFilesIn {}/*R1_001.fastq.gz {}/*R2_001.fastq.gz --runThreadN 6 --alignIntronMin 70 --alignIntronMax 100000 --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 1 --outFilterMultimapScoreRange 0 --outFilterMismatchNmax 5 --outFileNamePrefix {}/ --outStd Log --readFilesCommand zcat"
ls -d Sample_*[0-9] | xargs -I {} sh -c "mv {}/*bam {}.bam"
ls -d Sample_*[0-9] | parallel "samtools index {}.bam"
Rscript rnaseq_read_count_1.R
