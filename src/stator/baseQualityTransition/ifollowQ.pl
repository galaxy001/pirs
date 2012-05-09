#!/usr/bin/env perl
use strict;
use warnings;

die "Usage: $0 <LENtoStat> <read1.fq.gz> <outprefix>.fQ{out,dat}\n" if @ARGV != 3;
my ($LENtoStat,$fq1,$out)=@ARGV;
my $Qchr = 64;
if ($LENtoStat<0) {
	$Qchr = 33;
	$LENtoStat = -$LENtoStat;
}
$out .= ".$LENtoStat";

sub openfile($) {
    my ($filename)=@_;
    my $infile;
    if ($filename=~/.bz2$/) {
	    open( $infile,"-|","bzip2 -dc $filename") or die "Error opening $filename: $!\n";
    } elsif ($filename=~/.gz$/) {
     	open( $infile,"-|","gzip -dc $filename") or die "Error opening $filename: $!\n";
    } else {open( $infile,"<",$filename) or die "Error opening $filename: $!\n";}
    return $infile;
}
sub readfq($) {
	my $fh=$_[0];
	defined(my $a=<$fh>) or return [1];
	chomp($a);
	chomp(my $b=<$fh>) or return [2];
	chomp(my $c=<$fh>) or return [3];
	chomp(my $d=<$fh>) or return [4];
#print "$b\n$d\n";
	$b=substr($b,0,length($d)) if $d =~ s/B+$//;
#print "$a\n$b\n$c\n$d\n---\n";
	return [] unless length($d);
	return [$a,$b,$c,$d];
}
sub getQ($) {
    my @Qstr=split //,$_[0];
    my @Qvalue=();
    push @Qvalue,ord($_)-$Qchr for @Qstr;
    return \@Qvalue;
}
sub cal($) {
    my $Qhashes=$_[0];
    my ($x,$xx,$n,$cnt)=(0,0,0);
    my ($max,$min,$common,$maxcnt)=(0,(keys %$Qhashes)[0],0,0);
    for my $k (keys %$Qhashes) {
        $cnt = $$Qhashes{$k};
        $x += $k * $cnt;
        $xx += $k*$k * $cnt;
        $n += $cnt;
        $max = $k if $max < $k;
        $min = $k if $min > $k;
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

my ($ReadCnt,%statQ,$ret,%PlotReadsQavgHist)=(0);

#my $LENtoStat=10;
sub statQ($) {
    my $Qvalues=$_[0];
    my $Qlen=scalar @$Qvalues;
	my $Qsum=0;
	$Qsum += $_ for @$Qvalues;
	++$PlotReadsQavgHist{int(2*$Qsum/$Qlen)/2};
	++$PlotReadsQavgHist{-1};
    return "[x]Read Length must >= $LENtoStat." if $Qlen<$LENtoStat;
    for my $p (0..$Qlen-$LENtoStat-1) {
#   for my $p (19..$Qlen-$LENtoStat-1) {
        for my $q ($p+1..$p+$LENtoStat) {
            ++$statQ{$$Qvalues[$p]}{$$Qvalues[$q]};
        }
    }
}
my $FQa=openfile($fq1);
#my $FQb=openfile($fq2);

my $fqitem=&readfq($FQa);
while (@$fqitem != 1) {
    if (@$fqitem == 0) {
        $fqitem=&readfq($FQa);
        next;
    }
    my ($id,$Q)=@$fqitem[0,3];
    my $Qvalues=&getQ($Q);
    &statQ($Qvalues);
    $fqitem=&readfq($FQa);
    ++$ReadCnt;
}
close $FQa;

open OUT,'>',$out.'.fQout' or die "Error opening $out.fQout: $!\n";
open OD,'>',$out.'.fQdat' or die "Error opening $out.fQdat: $!\n";
print OUT "#ReadsCnt=$ReadCnt LENtoStat=$LENtoStat\n#",join("\t",qw/ Q cnt max min common mean std /),"\n";
print OD "#ReadsCnt=$ReadCnt LENtoStat=$LENtoStat\n#Q\toutMean\t",join("\t",(2..40)),"\n";
my ($above,$below,$at)=(0,0,0);
for my $k (sort {$a<=>$b} keys %statQ) {
    $ret=&cal($statQ{$k});
    print OUT join("\t",$k,@$ret),"\n";
    print OD "$k\t$$ret[-2]";
    for my $oq (2..40) {
        my $out='-';
        if (exists $statQ{$k}{$oq}){# && $k>2 && $oq>2) {
            $out=$statQ{$k}{$oq};
            if ($oq>$k) {$above+=$out;}# if $k>2;}
             elsif ($oq==$k) {$at+=$out;}
             else {$below+=$out;}
        }
        print OD "\t$out";
    }
    print OD "\n";
}
#print OUT "Y\t$_\t$Yrange{$_}\n" for sort {$a<=>$b} keys %Yrange;
close OUT;
print OD "# Above: $above\n# At: $at\n# Below: $below\n";

print OD "<<END\n\n";

my $MaxQ = 40;
print OD "[AvgQonReads]
#Total Reads Mean Quality values: $PlotReadsQavgHist{-1}
#Q\tCount\tRatio\n";
for my $q (0..2*$MaxQ) {
	$PlotReadsQavgHist{$q/2} = 0 unless exists $PlotReadsQavgHist{$q/2};
	print OD join("\t",$q/2,$PlotReadsQavgHist{$q/2},
		$PlotReadsQavgHist{$q/2}/$PlotReadsQavgHist{-1},"\n");
}
print OD "<<END\n";
close OD;

