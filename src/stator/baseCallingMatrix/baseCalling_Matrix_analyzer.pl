#!/bin/env perl

=head1 Name

	error_matrix_analyzer.pl

=head1 Description

	analyse the matrix about reads'' characteristic count from comparison results

	the program need the matrix file as input.It will creat 7 files:
	  *.misVSerr.base.stat : mismatch rate and error rate calculate by quality value for every cycle
	  *.qualVSmis.stat : compair the real mismatch rate and the calculated error rate for each quality value
	  *.transform.cycle.stat : the rate of the reference nucleotide be sequenced to each nucleotide for every cycle
	  *.transform.qual.stat : the rate of the reference nucleotide be sequenced to each nucleotide for every quality
	  *.transform.avg.stat : the average the reference nucleotide be sequenced to each nucleotide for all cycle along read
	  *.qual.mat.dis : the distribution of quality value for matched read base
	  *.qual.mis.dis : the distribution of quality value for mismatched read base
	  *.err2mis : from released/pIRS/bwasam/matrixsummer.pl

=head1 Version
	
	Author: Shi Yujian <shiyujian@genomics.org.cn>, Hu Xuesong <galaxy001@gmail.com>
	Version: 1.21 , Date:20120302

=head1 Usage

	perl error_matrix_analyzer.pl [option]
	  -i <str>       matrix file
	  -o <str>       output prefix
	  -m <num>       min quality score[default:0]
	  -x <num>       max quality score[default:40]
	  -B             ignore the bases that quality is B or #
	  -h             help

=cut

use strict;
use warnings;
use Getopt::Long;

##get options from command line into variables and set default values
my ($len, $in_file, $out, $min_qual, $max_qual, $help, $maskB);
GetOptions(
	"i:s"=>\$in_file,
	"o:s"=>\$out,
	"m:i"=>\$min_qual,
	"x:i"=>\$max_qual,
	"B"=>\$maskB,
	"help"=>\$help
);

die `pod2text $0` if($help || !defined $in_file);
$min_qual ||= 0;
$max_qual ||= 40;
if($in_file=~/\/([^\/]+)$/)
{
	$out ||= $1;
}
else
{
	$out ||= $in_file;
}

my %seq2bit = ("A"=>0,"C"=>1,"G"=>2,"T"=>3,"_All"=>4);
my @bit2seq = ("A","C","G","T","_All");
my @matrix;
open IN,$in_file || die"$!";

while(<IN>)
{
	if(/\[DistMatrix\]/)
	{
		$_=<IN>;
		die "input file format error" unless(/#Ref/);
		last;
	}
}

while(<IN>)
{
	next if(/#/);
	next unless(/\S+/);
	last if(/END/);
	chomp;
	my @line = split;
	my $k = 2;
	for(my $i=0;$i<4;$i++)
	{
		for(my $j=$min_qual;$j<=$max_qual;$j++)
		{
			$matrix[$seq2bit{$line[0]}][$line[1]-1][$i][$j] = $line[$k];
			$k++;
		}
		if($maskB)
		{
			$matrix[$seq2bit{$line[0]}][$line[1]-1][$i][2] = 0;
		}
	}
}
close IN;

my @q2e;
for(my $i=$min_qual;$i<=$max_qual;$i++)
{
	$q2e[$i] = 10 ** (-$i / 10);
}

my(@sum_row,@sum_cycle,@sum_ref,@sum_cycle_match,@sum_cycle_err,@transform_cycle,@transform_avg,@transform_qual,@transform_qual_sum);
my(@match_qual_cycle,@mis_qual_cycle);
my(@match_qual_sum,$match_sum,@mis_qual_sum,$mis_sum,$err_sum);
my ($value,%MismatchBYQ,%MismatchBYQgenome);
for(my $ref=0;$ref<@matrix;$ref++)
{
	for(my $cycle=0;$cycle<@{$matrix[$ref]};$cycle++)
	{
		for(my $base=0;$base<@{$matrix[$ref][$cycle]};$base++)
		{
			my $sum_tmp = 0;
			for(my $qual=$min_qual;$qual<=$max_qual;$qual++)
			{
				$value = $matrix[$ref][$cycle][$base][$qual];
				$sum_tmp += $value;
				my $err_tmp = $value * $q2e[$qual];
				$sum_cycle_err[$cycle] += $err_tmp;
				$err_sum += $err_tmp;
				$transform_qual[$qual][$ref][$base] += $value;
				if($ref == $base)
				{
					$sum_cycle_match[$cycle] += $value;
					$match_qual_cycle[$cycle][$qual] += $value;
					$match_qual_sum[$qual] += $value;
					$match_sum += $value;
					$MismatchBYQ{$base}{$qual}->[1] +=$value;
					$MismatchBYQgenome{$ref}{$qual}->[1] +=$value;
					$MismatchBYQ{$seq2bit{'_All'}}{$qual}->[1] +=$value;
					$MismatchBYQgenome{$seq2bit{'_All'}}{$qual}->[1] +=$value;
				}
				else
				{
					$mis_qual_cycle[$cycle][$qual] += $value;
					$mis_qual_sum[$qual] += $value;
					$mis_sum += $value;
					$transform_qual_sum[$qual][$ref] += $value;
					$MismatchBYQ{$base}{$qual}->[0] +=$value;
					$MismatchBYQgenome{$ref}{$qual}->[0] +=$value;
					$MismatchBYQ{$seq2bit{'_All'}}{$qual}->[0] +=$value;
					$MismatchBYQgenome{$seq2bit{'_All'}}{$qual}->[0] +=$value;
				}
			}
			$sum_row[$ref][$cycle] += $sum_tmp;
			$transform_cycle[$cycle][$ref][$base] += $sum_tmp;
			$transform_avg[$ref][$base] += $sum_tmp;
		}
		$sum_cycle[$cycle] += $sum_row[$ref][$cycle];
		$sum_ref[$ref] += $sum_row[$ref][$cycle];
	}
}

#open OUT,">$out\.ratio" || die "$!";
my $total_base = $match_sum + $mis_sum;
#print OUT "#total_base_number:\t$total_base\n";
#print OUT "#base_rario:\t";
#for(my $ref=0;$ref<@matrix;$ref++)
#{
#	my $sum_ref_rate = $sum_ref[$ref] / $total_base;
#	print OUT "$bit2seq[$ref]:$sum_ref_rate\t";
#}
#print OUT "\n#ref_base\tcycle\t";
#for(my $base=0;$base<@{$matrix[0][0]};$base++)
#{
#	for(my $qual=$min_qual;$qual<=$max_qual;$qual++)
#	{
#		print OUT "$bit2seq[$base]$qual\t";
#	}
#}
#print OUT "\n";
#for(my $ref=0;$ref<@matrix;$ref++)
#{
#	for(my $cycle=0;$cycle<@{$matrix[$ref]};$cycle++)
#	{
#		my $cycle_tmp = $cycle + 1;
#		print OUT "$bit2seq[$ref]\t$cycle_tmp\t";
#		for(my $base=0;$base<@{$matrix[$ref][$cycle]};$base++)
#		{
#			for(my $qual=$min_qual;$qual<=$max_qual;$qual++)
#			{
#				my $ratio = 0;
#				$ratio = $matrix[$ref][$cycle][$base][$qual] / $sum_row[$ref][$cycle] if($sum_row[$ref][$cycle]);
#				print OUT "$ratio\t";
#			}
#		}
#		print OUT "\n";
#	}
#}
#close OUT;

open OUT,">$out\.misVSerr\.base\.stat" || die "$!";
my $avg_mis_rate = $mis_sum / $total_base;
my $avg_err_rate = $err_sum / $total_base;
print OUT "#total_base:\t$total_base\n";
print OUT "#avg_mis_rate:\t$avg_mis_rate\n";
print OUT "#avg_err_rate:\t$avg_err_rate\n";
print OUT "#cycle\tmis_rate\terr_rate\n";
for(my $cycle=0;$cycle<@sum_cycle;$cycle++)
{
	my $mis_rate_tmp = 1 - $sum_cycle_match[$cycle] / $sum_cycle[$cycle];
	my $err_rate_tmp = $sum_cycle_err[$cycle] / $sum_cycle[$cycle];
	my $cycle_tmp = $cycle + 1;
	print OUT "$cycle_tmp\t$mis_rate_tmp\t$err_rate_tmp\n";
}
close OUT;

open OUT,">$out\.qualVSmis\.stat" || die "$!";
print OUT "qual\tQ_err\tmis_rate\n";
for(my $qual=$min_qual;$qual<=$max_qual;$qual++)
{
	my $mis_tmp = 0;
	$mis_tmp = $mis_qual_sum[$qual] / ($match_qual_sum[$qual]+$mis_qual_sum[$qual]) if($match_qual_sum[$qual]+$mis_qual_sum[$qual]);
#	my $qual_tmp = chr($qual + 64);
	my $qual_tmp = $qual;
	print OUT "$qual_tmp\t$q2e[$qual]\t$mis_tmp\n";
}
close OUT;

open OUT,">$out\.transform\.cycle\.stat" || die "$!";
open OUT2,">$out\.transform\.avg\.stat" || die "$!";
print OUT "#cycle\t";
print OUT2 "ref\\read\t";
my $flag = 0;
for(my $ref=0;$ref<@transform_avg;$ref++)
{
	for(my $base=0;$base<@{$transform_avg[$ref]};$base++)
	{
#		next if($base == $ref);
		print OUT "$bit2seq[$ref]\-\>$bit2seq[$base]\t";
		if(!$flag)
		{
			print OUT2 "$bit2seq[$base]\t";
		}
	}
	if(!$flag)
	{
		print OUT2 "\n";
		$flag=1;
	}
}
print OUT "\n#avg\t";
for(my $ref=0;$ref<@transform_avg;$ref++)
{
	print OUT2 "$bit2seq[$ref]\t";
	for(my $base=0;$base<@{$transform_avg[$ref]};$base++)
	{
#		next if($base == $ref);
		my $transform_rate_tmp = 0;
#		$transform_rate_tmp = $transform_avg[$ref][$base] / ($sum_ref[$ref] - $transform_avg[$ref][$ref]) if(($sum_ref[$ref] - $transform_avg[$ref][$ref]));
		$transform_rate_tmp = $transform_avg[$ref][$base] / $sum_ref[$ref] if($sum_ref[$ref]);
		print OUT "$transform_rate_tmp\t";
		print OUT2 "$transform_rate_tmp\t";
	}
	print OUT2 "\n";
}
close OUT2;
print OUT "\n";
for(my $cycle=0;$cycle<@transform_cycle;$cycle++)
{
	my $cycle_tmp = $cycle + 1;
	print OUT "$cycle_tmp\t";
	for(my $ref=0;$ref<@transform_avg;$ref++)
	{
		for(my $base=0;$base<@{$transform_avg[$ref]};$base++)
		{
#			next if($base == $ref);
			my $transform_rate_tmp = 0;
#			$transform_rate_tmp = $transform_cycle[$cycle][$ref][$base] / ($sum_row[$ref][$cycle] - $transform_cycle[$cycle][$ref][$ref]) if(($sum_row[$ref][$cycle] - $transform_cycle[$cycle][$ref][$ref]));
			$transform_rate_tmp = $transform_cycle[$cycle][$ref][$base] / $sum_row[$ref][$cycle] if($sum_row[$ref][$cycle]);
			print OUT "$transform_rate_tmp\t";
		}
	}
	print OUT "\n";
}
close OUT;

open OUT,">$out\.transform.qual.stat" || die "$!";
print OUT "#qual\t";
for(my $ref=0;$ref<@transform_avg;$ref++)
{
	for(my $base=0;$base<@{$transform_avg[$ref]};$base++)
	{
		next if($ref == $base);
		print OUT "$bit2seq[$ref]\-\>$bit2seq[$base]\t";
	}
}
print OUT "\n#avg\t";
for(my $ref=0;$ref<@transform_avg;$ref++)
{
	for(my $base=0;$base<@{$transform_avg[$ref]};$base++)
	{
		next if($ref == $base);
		my $transform_rate_tmp = 0;
		$transform_rate_tmp = $transform_avg[$ref][$base] / ($sum_ref[$ref] - $transform_avg[$ref][$ref]) if($sum_ref[$ref] - $transform_avg[$ref][$ref]);
#		$transform_rate_tmp = $transform_avg[$ref][$base] / $sum_ref[$ref] if($sum_ref[$ref]);
		print OUT "$transform_rate_tmp\t";
	}
}
print OUT "\n";
for(my $qual=$min_qual;$qual<=$max_qual;$qual++)
{
	print OUT "$qual\t";
	for(my $ref=0;$ref<@transform_avg;$ref++)
	{
		for(my $base=0;$base<@{$transform_avg[$ref]};$base++)
		{
			next if($ref == $base);
			my $transform_rate_tmp = 0;
			$transform_rate_tmp = $transform_qual[$qual][$ref][$base] / $transform_qual_sum[$qual][$ref] if($transform_qual_sum[$qual][$ref]);
#			$transform_rate_tmp = $transform_qual[$qual][$ref][$base] / ($match_qual_sum[$qual]+$mis_qual_sum[$qual]) if($match_qual_sum[$qual]+$mis_qual_sum[$qual]);
			print OUT "$transform_rate_tmp\t";
		}
	}
	print OUT "\n";
}
close OUT;

open OUT1,">$out\.qual\.mat\.dis" || die "$!";
open OUT2,">$out\.qual\.mis\.dis" || die "$!";
print OUT1 "#qual\tavg\t";
print OUT2 "#qual\tavg\t";
for(my $cycle=0;$cycle<@match_qual_cycle;$cycle++)
{
	my $cycle_tmp = $cycle + 1;
	print OUT1 "cycl:$cycle_tmp\t";
	print OUT2 "cycl:$cycle_tmp\t";
}
print OUT1 "\n";
print OUT2 "\n";
for(my $qual=$min_qual;$qual<=$max_qual;$qual++)
{
#	my $qual_tmp = chr($qual+64);
	my $match_tmp = $match_qual_sum[$qual] / $match_sum;
	my $mis_tmp = $mis_qual_sum[$qual] / $mis_sum;
	print OUT1 "$qual\t$match_tmp\t";
	print OUT2 "$qual\t$mis_tmp\t";
	for(my $cycle=0;$cycle<@match_qual_cycle;$cycle++)
	{
		$match_tmp = 0;
		$match_tmp = $match_qual_cycle[$cycle][$qual] / $sum_cycle_match[$cycle] if($sum_cycle_match[$cycle]);
		my $mis_tmp = 0;
		$mis_tmp = $mis_qual_cycle[$cycle][$qual] / ($sum_cycle[$cycle] - $sum_cycle_match[$cycle]) if($sum_cycle[$cycle] - $sum_cycle_match[$cycle]);
		print OUT1 "$match_tmp\t";
		print OUT2 "$mis_tmp\t";
	}
	print OUT1 "\n";
	print OUT2 "\n";
}
close OUT1;
close OUT2;

sub misRorTwo($$) {
	my ($mismatch,$match)=@_;
	if ($mismatch+$match==0) {
		if ($mismatch == 0) {
			return 0;
		} else {
			return 2;
		}
	} else {
		return $mismatch/($mismatch+$match);
	}
}
#ddx \%MismatchBYQ;
open OA,'>',$out.'.err2mis' or die "Error: $!\n";
print OA "#Read\tQ\tErrRate\tMismatchRate\tbyRefMismatchRate\tmismatch,match\n";
for my $read (sort keys %MismatchBYQ) {
    for my $Q (sort {$a<=>$b} keys %{$MismatchBYQ{$read}}) {
        my ($mismatch,$match)=@{$MismatchBYQ{$read}{$Q}};
		my ($gmismatch,$gmatch)=@{$MismatchBYQgenome{$read}{$Q}};
        next unless ($match+$gmatch+$mismatch+$gmismatch);
        print OA "$bit2seq[$read]\t$Q\t",10**(-$Q/10),"\t",misRorTwo($mismatch,$match),"\t",misRorTwo($gmismatch,$gmatch),"\t$mismatch,$match\t$gmismatch,$gmatch","\n";
    }
}
close OA;
