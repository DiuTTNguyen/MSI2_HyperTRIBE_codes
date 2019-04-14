This is a data processing pipeline for HyperTRIBE paired-end RNA-Seq data. It
is mostly based on the Broad pipeline for calling variants in RNA-Seq
(https://software.broadinstitute.org/gatk/documentation/article.php?id=3891)
with minor changes.

Assuming you get your fastq files from IGO, files from each sample should be
in a separate directory.  Please choose the right organism and then run 
tribe_alignment.sh, which will produce the BAM alignments. After that, run
variant_calling.sh and you will get VCF files containing mutations from each
library. Finally, you can run 'Rscript rnaseq_read_count.R' and it will produce
an Excel worksheet containing the full and annotated allele read count table.
Note that each step in variant calling will produce a separate bam
file. To avoid running out of space, you should remove the temporary files
that you no longer need.

Software used:
STAR (2.6.0a)
java (1.8.0_131)
picard (2.7.11)
GATK (3.8-1-0)
R (3.5.0)

Genome index and dbSNP annotation on lilac:
ln -s /data/leslie/luyuheng/index/GRCh38_Gencode25 hg19_star
ln -s /data/leslie/luyuheng/index/GRCm38_Gencode25 mm10_star
ln -s /home/chinc1/hg19 hg19_anno
ln -s /data/leslie/luyuheng/index/dbsnp_mouse dbsnp_mouse
