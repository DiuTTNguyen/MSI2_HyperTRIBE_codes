#BSUB-J variant_call        # name of the job
#BSUB-W 6:00                   	# Time limit 
#BSUB-n 9
#BSUB-R "span[ptile=9]"
#BSUB-M 30

#BSUB-o stdout
#BSUB-o stderr

# Load Java 8
module add java

# Indexing the fasta sequence
# java -jar ~/jartools/picard.jar CreateSequenceDictionary REFERENCE=mm10_star/GRCm38.primary_assembly.genome.fa OUTPUT=mm10_star/GRCm38.primary_assembly.genome.dict 
# samtools faidx mm10_star/GRCm38.primary_assembly.genome.fa
# Sort dbSNP VCF according to reference sequence
# java -jar ~/jartools/picard.jar SortVcf I=dbsnp_mouse/dbsnp_151_mouse.vcf  O=dbsnp_mouse/dbsnp_151_mouse.sorted.vcf  SEQUENCE_DICTIONARY=mm10_star/GRCm38.primary_assembly.genome.dict 

# Picard: add RG & mark duplicates
ls -d Sample_*[0-9] | parallel 'java -jar ~/jartools/picard.jar AddOrReplaceReadGroups I={}.bam O={}.rg_added.bam RGID=id RGLB=library RGPL=ILLUMINA RGPU=machine RGSM=sample'
ls -d Sample_*[0-9] | parallel 'java -jar ~/jartools/picard.jar MarkDuplicates I={}.rg_added.bam O={}.dedupped.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics'
# GATK: Split reads
ls -d Sample_*[0-9] | parallel 'java -jar ~/jartools/GenomeAnalysisTK.jar -T SplitNCigarReads -R mm10_star/GRCm38.primary_assembly.genome.fa -I {}.dedupped.bam -o {}.split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS'
# GATK: Call SNPs & run QC
ls -d Sample_*[0-9] | parallel 'java -jar ~/jartools/GenomeAnalysisTK.jar -T HaplotypeCaller -R mm10_star/GRCm38.primary_assembly.genome.fa -I {}.split.bam -dontUseSoftClippedBases -stand_call_conf 10.0 -o {}.HaplotypesR.vcf --dbsnp dbsnp_mouse/dbsnp_151_mouse.sorted.vcf'
ls -d Sample_*[0-9] | parallel 'java -jar ~/jartools/GenomeAnalysisTK.jar -T VariantFiltration -R mm10_star/GRCm38.primary_assembly.genome.fa -V {}.HaplotypesR.vcf -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o {}.FinalR.vcf'
