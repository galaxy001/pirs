#!/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);

#main:
open I,'<','ifQmerge.dat' or die "Error opening ifQmerge.dat: $!\n";
open OUT,'>','hist.dat' or die "Error opening hist.dat: $!\n";
while(<I>) {
	if (/^\[AvgQonReads\]\n/) {
		while(<I>){
			print OUT $_ unless /^<<END/;
		}
	}
}
close I;
close OUT;

system('gnuplot',"$Bin/hist.dem");
