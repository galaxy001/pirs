#!/usr/bin/env perl
use strict;
use warnings;

die "Usage: $0 <read1.fq.gz> <out>\n" if @ARGV != 2;
my ($fq1,$out)=@ARGV;

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
	defined(my $a=<$fh>) or return [];
	chomp($a);
	chomp(my $b=<$fh>) or return [];
	chomp(my $c=<$fh>) or return [];
	chomp(my $d=<$fh>) or return [];
	return [$a,$b,$c,$d];
}
sub getQ($) {
    my @Qstr=split //,$_[0];
    my @Qvalue=();
    push @Qvalue,ord($_)-64 for @Qstr;
    return \@Qvalue;
}
sub statQ($) {
    my $Qvalues=$_[0];
    my $Qlen=scalar @$Qvalues;
    die "[x]Read Length must >= 50." if $Qlen<50;
    my ($Q1,$Q5,$Q10,$Q15,$Qs15,$Q30,$Q50,$Qall);
    $Q1=$$Qvalues[0];
    $Q5+=shift @$Qvalues for (1..5);
    #$Q10=$Q15=$Qs15=$Q30=$Q50=$Qall=$Q5;
    $Q10=$Q5;
    $Q10+=shift @$Qvalues for (1..5);
    $Q15=$Q10;
    $Q15+=shift @$Qvalues for (1..5);
    $Qs15+=shift @$Qvalues for (1..15);
    $Q30=$Q15+$Qs15;
    $Q50=$Q30;
    $Q50+=shift @$Qvalues for (1..20);
    $Qall=$Q50;
    $Qall+=shift @$Qvalues while (@$Qvalues);
    return [$Q1,$Q5/5,$Q10/10,$Q15/15,$Qs15/15,$Q30/30,$Q50/50,$Qall/$Qlen];
}
my @statQlabel=qw/Q1 Q5 Q10 Q15 Q16t30 Q30 Q50 Qall/;

my $FQa=openfile($fq1);
#my $FQb=openfile($fq2);

my (%ReadsinTile,%Xrange,%Yrange,%statQ);

my $fqitem=&readfq($FQa);
while (@$fqitem > 0) {
    my ($id,$Q)=@$fqitem[0,3];
    my ($Tile,$X,$Y)=(split /[:#]/,$id)[2,3,4];
    my $Qvalues=&getQ($Q);
    my $Qstat=&statQ($Qvalues);
    #my ($Surface,$Swath,$subTile)=(split //,$Tile)[0,1,3];
    #print "$Surface,$Swath,$subTile\t$Tile,$X,$Y\n";
    ++$ReadsinTile{$Tile};
    $statQ{$Tile}=[] unless exists $statQ{$Tile};
    my $i=0;
    for my $Qv (@$Qstat) {
        ${$statQ{$Tile}}[$i] += $Qv;
        ++$i;
    }
    ++$Xrange{$X}; ++$Yrange{$Y};
    $fqitem=&readfq($FQa);
}
close $FQa;

open OUT,'>',$out or die "Error opening $out: $!\n";
for my $tile (sort keys %ReadsinTile) {
    print OUT "Tile\t$tile\t$ReadsinTile{$tile}\n";
}
print OUT "#Qstat\tTile\t",join("\t",@statQlabel),"\n";
for my $tile (sort keys %statQ) {
    $_ /=$ReadsinTile{$tile} for @{$statQ{$tile}};
    print OUT "Qstat\t$tile\t",join("\t",@{$statQ{$tile}}),"\n";
}
print OUT "X\t$_\t$Xrange{$_}\n" for sort {$a<=>$b} keys %Xrange;
print OUT "Y\t$_\t$Yrange{$_}\n" for sort {$a<=>$b} keys %Yrange;
close OUT;



