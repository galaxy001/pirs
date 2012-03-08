#!/bin/env perl
use strict;
use warnings;

my $SAMTOOLSBIN="/ifs1/ST_ASMB/USER/yuanjy/huxuesong/tmp/group/rev/test/samtools";
my $MAXREADStoCHECK=10000;

die "Usage: $0 <single_bam_sam_file> <output>\n" if @ARGV<2;
my $name=shift;
my $out=shift;
my ($READLEN,$lines)=(0,0);
if ($name =~ /\.bam$/) {
	open IN,'-|',"$SAMTOOLSBIN view -f 3 -F 1536 $name" or die "Error opening $name : $!\n";
} elsif ($name =~ /\.sam\.gz$/) {
	open IN,'-|',"$SAMTOOLSBIN view -f 3 -F 1536 -S $name" or die "Error opening $name : $!\n";
} elsif ($name =~ /\.sam$/) {
	open IN,'-|',"$SAMTOOLSBIN view -f 3 -F 1536 -S $name" or die "Error opening $name : $!\n";
} else {
	die "[x]Unsupport file type.";
}

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
$lines=0;
if ($name =~ /\.bam$/) {
	open IN,'-|',"$SAMTOOLSBIN view -f 3 -F 1536 $name" or die "Error opening $name : $!\n";
} elsif ($name =~ /\.sam\.gz$/) {
	open IN,'-|',"$SAMTOOLSBIN view -f 3 -F 1536 -S $name" or die "Error opening $name : $!\n";
} elsif ($name =~ /\.sam$/) {
	open IN,'-|',"$SAMTOOLSBIN view -f 3 -F 1536 -S $name" or die "Error opening $name : $!\n";
} else {
	die "[x]Unsupport file type.";
}

my (%Cnt,%DistIns,%DistDel);
my ($Read12,$PosShift,$cigar,$DelShift);
while (<IN>) {
	next if /^@\w\w\t\w\w:/;
	chomp;
	my @read1=split /\t/;
	#next unless ($read1[1] & 3) == 3;  # paired + mapped in a proper pair
	#next if $read1[1] >= 256;   # not primary || QC failure || optical or PCR duplicate
	next if $read1[5] =~ /S/;	# only IDM, no S.
	my $OPT = join("\t",@read1[11 .. $#read1]);
	next if $OPT =~ /\bXT:A:R\b/;
	if ($read1[1] & 64) {
		$Read12 = 1;
	} elsif ($read1[1] & 128) {
		$Read12 = 2;
	} else {
		warn "[w]",join("\t",@read1),"\n";
		next;
	}
	++$Cnt{$Read12}{'All'};
	if ($read1[1] & 16) {	#  | r  | 0x0010 | strand of the query (1 for reverse)   |
		$PosShift = -$READLEN - 1;
		$DelShift = 0;
	} else {	# +
		$PosShift = 0;
		$DelShift = -1;
	}
	$cigar=$read1[5];
	if ($read1[5] =~ /(\d+)[ID]/) {	# I/D in reads
      # http://davetang.org/muse/2011/01/28/perl-and-sam/
      my $position = '1';
      while ($cigar !~ /^$/){
         if ($cigar =~ /^([0-9]+[MIDS])/){
            my $cigar_part = $1;
            if ($cigar_part =~ /(\d+)M/){
               $position += $1;
            } elsif ($cigar_part =~ /(\d+)I/){
				$Cnt{$Read12}{'Ins'} += $1;
				for my $p ($position .. ($position + $1 -1)) {
					$DistIns{$Read12}{abs($p+$PosShift)}++;
warn "$position -> ",$position+$PosShift,"\t$1\t$cigar\t$cigar_part\n$_\n" if abs($position+$PosShift)<=1 or abs($position+$PosShift)>=$READLEN;
				}
#warn "$position -> ",$position+$PosShift,"\t$cigar\t$cigar_part\n$_\n" if abs($position+$PosShift)<1 or abs($position+$PosShift)>$READLEN;
               $position += $1;
            } elsif ($cigar_part =~ /(\d+)D/){
				$Cnt{$Read12}{'Del'} += $1;
				my $p=abs($position+$DelShift+$PosShift);
				$DistDel{$Read12}{$p}{$1}++;	# 99M1D1M: D@99, not 100.
				$DistDel{$Read12}{-1}++;
warn "$position -> ",$p,"\t$1\t$cigar\t$cigar_part\n$_\n" if $p<1 or $p>=$READLEN;
               #$position += $1;
            } elsif ($cigar_part =~ /(\d+)S/){
               die "[!]Not ready for this!\n";
               #my $insertion = 'x' x $1;
               #substr($new_ref,$position,0,$insertion);
               #$position += $1;
            }
            $cigar =~ s/$cigar_part//;
         } else {
            die "Unexpected cigar: $cigar\n";
         }
      }
	}
	++$lines;
	last if $lines > 100*$MAXREADStoCHECK;
}
close IN;

open O,'>',$out or die "Error opening ${out} : $!\n";
print O "File=$name\nRead_Length=$READLEN\n\nRead\tType\tCount\tRatio\n";
for my $Read12 (sort keys %Cnt) {
	for (sort keys %{$Cnt{$Read12}}) {
		#next if $_ eq 'All';
		print O "$Read12\t$_\t$Cnt{$Read12}{$_}\t",$Cnt{$Read12}{$_}/$Cnt{$Read12}{'All'},"\n";
	}
}

print O "\nRead\tCycle\tDel\tCount\tRatio\n";
for my $Read12 (sort keys %DistDel) {
	for my $cyc (sort {$a<=>$b} keys %{$DistDel{$Read12}}) {
		next if $cyc == -1;
		for my $ins (sort {$a<=>$b} keys %{$DistDel{$Read12}{$cyc}}) {
			print O "$Read12\t$cyc\t$ins\t$DistDel{$Read12}{$cyc}{$ins}\t",$DistDel{$Read12}{$cyc}{$ins}/$DistDel{$Read12}{-1},"\n";
		}
	}
}

print O "\nRead\tCycle\tIns\tRatio\n";
for my $Read12 (sort keys %DistIns) {
	for my $cyc (sort {$a<=>$b} keys %{$DistIns{$Read12}}) {
		print O "$Read12\t$cyc\t$DistIns{$Read12}{$cyc}\t",$DistIns{$Read12}{$cyc}/$Cnt{$Read12}{'Ins'},"\n";
	}
}

close O;
