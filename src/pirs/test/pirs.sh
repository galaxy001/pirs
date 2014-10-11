#!/bin/bash
./../pirs diploid ref_seq.fa -s 0.001 -R 2 -d 0.0001 -v 0.000001 -c gzip -o ref22 >SimDiploid.out 2>SimDiploid.err

PROFILES="-B ../Profiles/Base-Calling_Profiles/humNew.PE100.matrix.gz -I ../Profiles/InDel_Profiles/phixv2.InDel.matrix -G ../Profiles/GC-depth_Profiles/humNew.gcdep_100.dat"

./../pirs simulate -A dist ref_seq.fa ref22.snp.indel.inversion.fa.gz -d -M lowercase -m 800 -l 100 -x 10 -v 40 -z ${PROFILES} -o Illumina >SimReads.out 2>SimReads.err
./../pirs simulate ref_seq.fa ref22.snp.indel.inversion.fa.gz -d -M lowercase -m 800 -l 100 -x 10 -v 40 -z ${PROFILES} -o EAMSS2 >SimReadsEAMSS2.out 2>SimReadsEAMSS2.err
