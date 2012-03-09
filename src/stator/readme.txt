illumia_reads_parameter_stator
	by Hu Xuesong @ BGI <huxuesong@genomics.org.cn>

GC%-Depth Profile Stat:
1. run soap and soap.coverage to get .depthsingle file(s). gzip is OK to over it.
2. run gc_coverage_bias on all depthsingle files. You will get gc-depth stat by 1 GC% and other files.
3. run gc_coverage_bias_plot on the gc-depth stat file. You'll get PNG plot.

Base-Calling Profile Stat:
1. run soap or bwa to get .{soap,single} or .sam/.bam file(s).
2. run error_matrix_calculator on those file(s). You will get *.count.matrix .
   [Sorted BAM/SAM files require more memory and thus not suggested.]
3. You can use error_matrix_merger to merge several .count.matrix files.
   However, it is up to you to keep the read length matches.

InDel Profile Stat:
1. choose samples with NO polymorphism InDel, such as the Coliphage samples that shipped with Illumina Sequencers.
2. run bwa to get .sam/.bam file.
3. run indelstat_sam_bam to get the profile.

Insert size & mapping ratio stat:
1. run soap or bwa to get .{soap,single} or .sam file(s).
2. run alignment_stator *.
* alignment_stator cannot stat. mapping ratio for sam files now.
[gzipped sam files can be read as STDIN with "zcat file.sam.gz|",
 bam files can be used as STDIN with "samtools view -f 3 -F 1792 -h file.bam|".]
