#!/bin/env perl

use strict;
use warnings;

die "Usage: perl $0 <matrix list> <merged outfile>" if(@ARGV<2);

open IN,$ARGV[0] || die "$!";
my(@handles,$HAND);
my $h=0;
while(<IN>)
{
	open $handles[$h++],$_ || die "$!";
}
close IN;

open OUT,">$ARGV[1]" || die "$!";

#merge the head of the input files
my $user="";
my(@map_reads,@map_base,@stat_base,@stat_reads,@mis_base);
my($map_reads_sum,$map_base_sum,$stat_base_sum,$stat_reads_sum,$read_len,$mis_base_sum)=(0,0,0,0,0,0);
my($ref,$base,$cycle,$qual)=(0,0,0,0);
my($qb_base,$qb_mis)=(0,0);
my($A,$C,$G,$T);
for(my $i=0;$i<@handles;$i++)
{
	$HAND = $handles[$i];
	my $data=<$HAND>;
	chomp $data;
	if($data=~/by\s+([\s\S]+)$/)
	{
		my $tmp=$1;
		$user .= ";" . $tmp if(!($user=~/$tmp/));
	}
	else{die "35\tplease check the input files carefully";}

	$data=<$HAND>;
	if($data=~/#Input\D+(\d+)\D+(\d+)/)
	{
		$map_reads[$i] = $1;
		$map_reads_sum += $1;
		$map_base[$i] = $2;
		$map_base_sum += $2;
	}
	else{die "45\tplease check the input files carefully";}

	$data=<$HAND>;
	if($data=~/#Total\D+(\d+)\D+(\d+)\D+(\d+)/)
	{
		die "read_len is not same" if($read_len !=0 && $read_len != $3);
		$stat_base[$i] = $1;
		$stat_base_sum += $1;
		$stat_reads[$i] = $2;
		$stat_reads_sum += $2;
		$read_len = $3;
	}
	else{die "57\tplease check the input files carefully";}

	$data=<$HAND>;
	if($data=~/#Dimensions\D+(\d+)\D+(\d+)\D+(\d+)\D+(\d+)/)
	{
		die "ref base classes not same" if($ref!=0 && $ref!=$1);
		die "cycle number not same" if($cycle!=0 && $cycle!=$2);
		die "seq base classes not same" if($base!=0 && $base!=$3);
		die "quality number not same" if($qual!=0 && $qual!=$4);
		$ref=$1;$cycle=$2;$base=$3;$qual=$4;
	}
	else{die "68\tplease check the input files carefully";}

	$data=<$HAND>;
	if($data=~/#Mismatch_base\D+(\d+)\D+([\d\.]+)/)
	{
		$mis_base[$i]=$1;
		$mis_base_sum += $1;
	}
	else{die "76\tplease check the input files carefully";}

	$data=<$HAND>;
	if($data=~/#QB_Bases\D+(\d+)\D+(\d+)/)
	{
		$qb_base += $1;
		$qb_mis += $2;
	}
	else{die "84\tplease check the input files carefully";}

	$data=<$HAND>;
	if($data=~/#Reference\D+([\d\.]+)\D+([\d\.]+)\D+([\d\.]+)\D+([\d\.]+)/)
	{
		$A += $1 * $stat_base[$i];
		$C += $2 * $stat_base[$i];
		$G += $3 * $stat_base[$i];
		$T += $4 * $stat_base[$i];
	}
	else{die "94\tplease check the input files carefully";}

	while(<$HAND>)
	{
		last if(/DistMatrix/);
	}
}

my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
my $date=sprintf "%02d:%02d:%02d,%4d-%02d-%02d",$hour,$min,$sec,$year+1900,$mon+1,$mday;
print OUT "#Generate \@ $date by$user\n";
print OUT "#Input [sam] file of mapped Reads: $map_reads_sum , mapped Bases $map_base_sum (no base stat for sam files)\n";
print OUT "#Total statistical Bases: $stat_base_sum , Reads: $stat_reads_sum of length $read_len\n";
print OUT "#Dimensions: Ref_base_number $ref, Cycle_number $cycle, Seq_base_number $base, Quality_number $qual\n";
my $mis_rate = $mis_base_sum / $stat_base_sum * 100;
print OUT "#Mismatch_base: $mis_base_sum, Mismatch_rate: $mis_rate %\n";
print OUT "#QB_Bases: $qb_base, QB_Mismatches: $qb_mis (bases with quality <= 2)\n";
$A /= $stat_base_sum;
$C /= $stat_base_sum;
$G /= $stat_base_sum;
$T /= $stat_base_sum;
print OUT "#Reference Base Ratio in reads: A $A %;   C $C %;   G $G %;   T $T %;\n";
print OUT "\n\[DistMatrix\]\n";

#my $type = "";
my $head = "";
for(my $i=0;$i<@handles;$i++)
{
	$HAND = $handles[$i];
	my $data = <$HAND>;
#	if($data=~/Type\s+\=\s+(\w+)/)
#	{
#		die "types are not same" if($type ne "" && $type ne $1);
#		$type = $1;
#	}
#	else{die "129\tplease check the input files carefully";}
	
#	$data = <$HAND>;
	die "the heads are not same" if($head ne "" && $head ne $data);
	$head = $data;
}
#print OUT "Type = $type\n$head";
print OUT $head;

$HAND = $handles[0];
while(<$HAND>)
{
	chomp;
	if(/END/)
	{
		for(my $i=1;$i<@handles;$i++)
		{
			$HAND = $handles[$i];
			my $data = <$HAND>;
			die "147\tplease check the input files carefully" unless($data=~/END/);
		}
		last;
	}
	my @sum = split;
	for(my $i=1;$i<@handles;$i++)
	{
		$HAND = $handles[$i];
		my $data = <$HAND>;
		my @line = split /\s+/,$data;
		die "157\tplease check the input files carefully" if($line[0] ne $sum[0] || $line[1] != $sum[1] || $#sum != $#line);
		for(my $j=2;$j<@sum;$j++)
		{
			$sum[$j] += $line[$j];
		}
	}
	for(my $j=0;$j<$#sum;$j++)
	{
		print OUT "$sum[$j]\t";
	}
	print OUT "$sum[-1]\n";
	$HAND = $handles[0];
}
print OUT "<<END\n";

for(my $i=0;$i<@handles;$i++)
{
	$HAND = $handles[$i];
	while(<$HAND>)
	{
		last if(/QTransMatrix/);
	}
}
print OUT "\n\[QTransMatrix\]\n";
#$type = "";
$head = "";
for(my $i=0;$i<@handles;$i++)
{
	$HAND = $handles[$i];
	my $data = <$HAND>;
#	if($data=~/Type\s+=\s+(\w+)/)
#	{
#		die "types are not same" if($type ne "" && $type ne $1);
#		$type = $1;
#	}
#	else{die "192\tplease check the input files carefully";}
	
#	$data = <$HAND>;
	die "the heads are not same" if($head ne "" && $head ne $data);
	$head = $data;
}
#print OUT "Type = $type\n$head";	
print OUT $head;

$HAND = $handles[0];
while(<$HAND>)
{
	chomp;
	if(/END/)
	{
		for(my $i=1;$i<@handles;$i++)
		{
			$HAND = $handles[$i];
			my $data = <$HAND>;
			die "210\tplease check the input files carefully" unless($data=~/END/);
		}
		last;
	}
	my @sum = split;
	for(my $i=1;$i<@handles;$i++)
	{
		$HAND = $handles[$i];
		my $data = <$HAND>;
		my @line = split /\s+/,$data;
		die "220\tplease check the input files carefully" if($line[0] ne $sum[0] || $line[1] != $sum[1] || $#sum != $#line);
		for(my $j=2;$j<@sum;$j++)
		{
			$sum[$j] += $line[$j];
		}
	}
	for(my $j=0;$j<$#sum;$j++)
	{
		print OUT "$sum[$j]\t";
	}
	print OUT "$sum[-1]\n";
	$HAND = $handles[0];
}
print OUT "<<END\n";

for(my $i=0;$i<@handles;$i++)
{
	$HAND = $handles[$i];
	close $HAND;
}
close OUT;
