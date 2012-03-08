./../pirs diploid -i ref_seq.fa -s 0.001 -a 2 -d 0.0001 -v 0.000001 -c 1 -o ref22 >SimDiploid.out 2>SimDiploid.err
./../pirs simulate -i ref_seq.fa -I ref22.snp.indel.invertion.fa.gz  -M 1 -m 800 -l 100 -x 10 -v 40 >SimReads.out 2>SimReads.err
