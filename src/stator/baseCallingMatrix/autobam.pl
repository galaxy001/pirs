#!/usr/bin/env perl
use strict;
use warnings;

my $SAMTOOLSBIN="/ifs1/ST_ASMB/USER/yuanjy/huxuesong/tmp/group/rev/test/samtools";
my $CALTUATORBIN="/ifs1/ST_ASMB/USER/yuanjy/huxuesong/tmp/group/rev/test/pl/baseCallingMatrix/baseCalling_Matrix_calculator.pl";
my $MAXREADStoCHECK=10000;

die "Usage: $0 <single_bamfile> <output_prefix> <other_options>\n" if @ARGV<3;
my $name=shift;
my $out=shift;
my ($READLEN,$lines)=(0,0);
open IN,'-|',"$SAMTOOLSBIN view -f 3 -F 1536 $name" or die "Error opening $name : $!\n";

while (<IN>) {
	next if /^@\w\w\t\w\w:/;
	chomp;
	my @read1=split /\t/;
	next unless ($read1[1] & 3) == 3;  # paired + mapped in a proper pair
	next if $read1[1] >= 256;   # not primary || QC failure || optical or PCR duplicate
	next unless $read1[5] =~ /^(\d+)M$/;
	$READLEN = $1 if $READLEN < $1;
	++$lines;
	last if $lines > $MAXREADStoCHECK;
}
close IN;

open O,'>',"${out}.log" or die "Error opening ${out}.log : $!\n";

my $CLI="$SAMTOOLSBIN view -f 3 -F 1792 $name | $CALTUATORBIN -bp sam -l $READLEN @ARGV -o $out >>${out}.log 2>${out}.err";
print "Read_Length of [$name]:\n$READLEN\n\nCLI:[$CLI]\n";
print O "BAM:[$name]\nRead_Length: $READLEN\n\nCLI:[$CLI]\n",'-' x 78,"\n\n";
system $CLI;

close O;
