#!/bin/env perl
use strict;
use warnings;
#use Data::Dump qw(dump ddx);

die "Usage: $0 <in.fQdat, ...>\n[!]Output filename is CONST.\n" if @ARGV < 1;
my ($MinQ,$MaxQ)=(2,40);

my ($ReadsCnt,$HistSum,%LentoStat,%Dat,%Hist)=(0,0);
for ($MinQ .. $MaxQ) {
	$Dat{$_}=[0];
}

sub readmatrix() {
	while(<>){
		next if /^#/;
		next if /^Type/;
		last if /<<END/;
		chomp;
		#s/\t\-(\t|$)/\t0\1/g;
		s/-/0/g;
		my @a=split /\t/;
		shift @a;	# base
		shift @a;	# cycle
		my $Qin=shift @a;
		my $Qlen=scalar @a;
		shift @a;	# 0
		shift @a;	# 1
		pop @a;	# RowSum
		#$MaxQ=$Qlen if $Qlen>$MaxQ;
		$Dat{$Qin}->[0] += $_ for @a;
		for my $q ($MinQ .. $MaxQ) {
			$Dat{$Qin}->[$q-$MinQ+1] += $a[$q-$MinQ];
		}
	}
	return;
}

sub cal($) {
    my @Qarray=@{$_[0]};
	my $n=shift @Qarray;
    my ($x,$xx,$cnt)=(0,0,0);
    my ($max,$min,$common,$maxcnt)=(0,$Qarray[0],0,0);
    for my $p (0..$#Qarray) {
		my $k=$p+$MinQ;
        $cnt = $Qarray[$p];
        $x += $k * $cnt;
        $xx += $k*$k * $cnt;
        #$n += $cnt;
		if ($cnt>0) {
			$max = $k if $max < $k;
			$min = $k if $min > $k;
		}
        if ($maxcnt<$cnt) {
            $maxcnt = $cnt;
            $common = $k;
        }
    }
    if ($n<2) {
        return [$n,0,0] if $n<1;
        return [$n,$max,0] if $n==1;
    }
    my $mean=$x/$n;
    my $std=($xx-$x*$mean)/$n-1;
    if ($std>0) {
        $std=sqrt($std);
    } elsif ($std==-1) {
        $std='i';
    } else {
        $std="sqrt($std)";
    }
    return [$n,$max,$min,$common,$mean,$std];
}

#main:
while(<>) {
	if (/^\[QTrans\]\n/) {
		readmatrix();
	}
}

my $LENtoStat = join(',',sort keys %LentoStat);
open OUT,'>','ifQmerge.out' or die "Error opening ifQmerge.out: $!\n";
open OD,'>','ifQmerge.dat' or die "Error opening ifQmerge.dat: $!\n";
print OUT "#ReadsCnt=$ReadsCnt LENtoStat=$LENtoStat\n#",join("\t",qw/ Q cnt max min common mean std /),"\n";
print OD "#ReadsCnt=$ReadsCnt LENtoStat=$LENtoStat\n#Q\toutMean\t",join("\t",(2..40)),"\n";
my ($above,$below,$at)=(0,0,0);
for my $k (sort {$a<=>$b} keys %Dat) {
    my $ret=&cal($Dat{$k});
	if ($$ret[0]) {
		print OUT join("\t",$k,@$ret),"\n";
		print OD "$k\t$$ret[-2]";
		for my $p (1..$#{$Dat{$k}}) {
			my $oq=$p-1+$MinQ;
			my $out='-';
			if (defined $Dat{$k}->[$p]){# && $k>2 && $oq>2) {
				$out=$Dat{$k}->[$p];
				if ($oq>$k) {$above+=$out;}# if $k>2;}
				 elsif ($oq==$k) {$at+=$out;}
				 else {$below+=$out;}
			}
			print OD "\t$out";
		}
		print OD "\n";
	}
}
#print OUT "Y\t$_\t$Yrange{$_}\n" for sort {$a<=>$b} keys %Yrange;
close OUT;
print OD "# Above: $above\n# At: $at\n# Below: $below\n";
print OD "<<END\n";
close OD;

for my $q (keys %Dat) {
	if ($Dat{$q}->[0] > 0) {
		$Dat{$q}->[$_] /= $Dat{$q}->[0] for (1 .. $MaxQ-$MinQ+1)
	} else {
		delete $Dat{$q};	#->[$_] = 0 for (1 .. $MaxQ-$MinQ+1)	
	}
}

open OD,'>','ifQplot.dat' or die "Error opening ifQplot.dat: $!\n";
for my $i ($MinQ..$MaxQ) {
	if (exists $Dat{$i}) {
		for my $j ($MinQ..$MaxQ) {
			print OD "$i\t$j\t$Dat{$i}->[$j-($MinQ-1)]\n"
		}
	} else {
		for my $j ($MinQ..$MaxQ) {
			print OD "$i\t$j\t0\n"
		}
	}
}
close OD;
