#!/usr/bin/perl -w

use strict;

die "Usage: perl $0 <gcdepth.txt> <span> <fit>" if(@ARGV<3);

open IN,$ARGV[0] || die "$!";
my(@x,@y,$span);
$span = $ARGV[1];
while(<IN>)
{
	chomp;
	next if(/^#/);
	my @line = split;
	push @x,$line[0];
	push @y,$line[1];
}
close IN;

my @f = @y;

for(my $i=0;$i<@x;$i++)
{
				my($a00,$a01,$a11,$d0,$d1)=(0,0,0,0,0);
				for(my $j=$i-$span;$j<=$i+$span;$j++)
				{
								if($j < 0 || $j >= @x){next;}
								my $w =  (1-(abs($x[$j]-$x[$i])/$span)**3)**1.5;
								$a00 += $w;
								$a01 += $w*$x[$j];
								$a11 += $w*$x[$j]*$x[$j];
								$d0 += $w*$y[$j];
								$d1 += $w*$x[$j]*$y[$j];
				}
				$f[$i] = ($a11*$d0-$a01*$d1)/($a11*$a00-$a01*$a01) + $x[$i] * ($a01*$d0-$a00*$d1)/($a01*$a01-$a00*$a11);
				if($f[$i] < 0){$f[$i] = 0;}
}


#for(my $i=$span;$i<@x-$span;$i++)
#{
#	my($a00,$a01,$a11,$d0,$d1)=(0,0,0,0,0);
#	for(my $j=$i-$span;$j<=$i+$span;$j++)
#	{
#		my $w =  (1-(abs($x[$j]-$x[$i])/$span)**3)**1.5;
#		$a00 += $w;
#		$a01 += $w*$x[$j];
#		$a11 += $w*$x[$j]*$x[$j];
#		$d0 += $w*$y[$j];
#		$d1 += $w*$x[$j]*$y[$j];
#	}
##	my $a0 = ($a11*$d0-$a01*$d1)/($a11*$a00-$a01*$a01);
##	my $a1 = ($a01*$d0-$a00*$d1)/($a01*$a01-$a00*$a11);
##	print "$a0\n$a1\n";
#	$f[$i] = ($a11*$d0-$a01*$d1)/($a11*$a00-$a01*$a01) + $x[$i] * ($a01*$d0-$a00*$d1)/($a01*$a01-$a00*$a11);
#}

open OUT,">$ARGV[2]" || die "$!";
print OUT "#x\ty\tf\n";
for(my $i=0;$i<@x;$i++)
{
	print OUT "$x[$i]\t$y[$i]\t$f[$i]\n";
}
close OUT;
