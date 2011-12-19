illumia_reads_parameter_stator
	by Hu Xuesong @ BGI <huxuesong@genomics.org.cn>

GC-Depth Stat:
1. run soap and soap.coverage to get .depthsingle file(s). gzip is OK to over it.
2. run gc_coverage_bias on all depthsingle files. You will get gc-depth stat by 1 GC% and other files.
3. run gc_coverage_bias_plot on the gc-depth stat file. You'll get PNG plot and a .gc file by 5 GC%.
4. Manually check the .gc file for any abnormal levels due to the lower depth on certain GC% windows.

Error Matrix Stat:
1. run soap or bwa to get .{soap,single} or .sam file(s).
2. run error_matrix_calculator on those file(s). You will get *.{count,ratio}.matrix .
3. You can use error_matrix_merger to merge several .{count,ratio}.matrix files.
   However, it is up to you to keep the read length matches.

Insert size & mapping ratio stat:
1. run soap or bwa to get .{soap,single} or .sam file(s).
2. run alignment_stator *.
* alignment_stator cannot stat. mapping ratio for sam files now.

