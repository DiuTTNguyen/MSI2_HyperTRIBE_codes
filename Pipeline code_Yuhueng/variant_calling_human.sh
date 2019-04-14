#BSUB-J MOLM13_variant_call        # name of the job
#BSUB-W 6:00                   	# Time limit 
#BSUB-n 9
#BSUB-R "span[ptile=9]"
#BSUB-M 30

#BSUB-o stdout
#BSUB-o stderr

# Load Java 8
module add java
# Picard: add RG & mark duplicates
ls -d Sample_*[0-9] | parallel 'java -jar ~/jartools/picard.jar AddOrReplaceReadGroups I={}.bam O={}.rg_added.bam RGID=id RGLB=library RGPL=ILLUMINA RGPU=machine RGSM=sample'
# Note: this step takes ~30G memory for each file; Be sure to reserve enough RAM
ls -d Sample_*[0-9] | xargs -I {} sh -c 'java -jar ~/jartools/picard.jar MarkDuplicates I={}.rg_added.bam O={}.dedupped.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics'
# hg19 genome we used for alignment is slightly different than the one used by picard 
# For now we have to rename and reorder the chrs to make them consistent
# Someone should try to fix this later
ls -d Sample_*[0-9] | parallel 'samtools reheader header.txt {}.dedupped.bam > {}.reheader.bam'
ls -d Sample_*[0-9] | parallel 'java -jar ~/jartools/picard.jar ReorderSam I={}.reheader.bam O={}.reordered.bam R=hg19_anno/ucsc.hg19.fasta CREATE_INDEX=TRUE'
# Split reads
ls -d Sample_*[0-9] | parallel 'java -jar ~/jartools/GenomeAnalysisTK.jar -T SplitNCigarReads -R hg19_anno/ucsc.hg19.fasta -I {}.reordered.bam -o {}.split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS'
# Call SNPs & annotate
ls -d Sample_*[0-9] | parallel 'java -jar ~/jartools/GenomeAnalysisTK.jar -T HaplotypeCaller -R hg19_anno/ucsc.hg19.fasta -I {}.split.bam -dontUseSoftClippedBases -stand_call_conf 10.0 -o {}.HaplotypesR.vcf --dbsnp hg19_anno/dbsnp_138.hg19.vcf'
ls -d Sample_*[0-9] | parallel 'java -jar ~/jartools/GenomeAnalysisTK.jar -T VariantFiltration -R hg19_anno/ucsc.hg19.fasta -V {}.HaplotypesR.vcf -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o {}.FinalR.vcf'
# Rscript rnaseq_read_count_1.R
