#!/usr/bin/env perl
use strict;
use warnings;

die "Usage: $0 <read1.fq.gz> <outprefix>.fqQ\n" if @ARGV != 2;
my ($fq1,$out)=@ARGV;
my $Qchr = 33;

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

my $QLen=0;

sub getQ($) {
    my @Qstr=split //,$_[0];
    my @Qvalue=();
    push @Qvalue,ord($_)-$Qchr for @Qstr;
	$QLen = $#Qvalue if $QLen < $#Qvalue;
    return \@Qvalue;
}

my ($ReadCnt,@statQv)=(0);

#my $LENtoStat=10;
sub statQ($) {
    my $Qvalues=$_[0];
    my $Qlen=$#$Qvalues;
	my $Qsum=0;
	$Qsum += $_ for @$Qvalues;
    for my $p (0..$Qlen) {
		$statQv[$p] += $$Qvalues[$p];
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

open OUT,'>',$out.'.fqQ' or die "Error opening $out.fQout: $!\n";
print OUT "# ReadsCnt: $ReadCnt\n# Cycle\tAvg_Q\tSum_Q\n";
for my $p (0..$QLen) {
	print OUT join("\t",$p+1,$statQv[$p]/$ReadCnt,$statQv[$p]),"\n";
}
close OUT;
