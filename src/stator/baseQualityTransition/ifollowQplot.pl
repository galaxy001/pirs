#!/usr/bin/env perl
use strict;
use warnings;

die "Usage: $0 <in.fQdat> > Output_to_STDOUT\n" if @ARGV != 1;
my ($in)=@ARGV;

my $MaxQ=40;

my %Dat;
open I,'<',$in or die $!;
while(<I>){
    next if /^#/;
    chomp;
    #s/\t\-(\t|$)/\t0\1/g;
    s/-/0/g;
    my @a=split /\t/;
    my $Qin=shift @a;
    my $Qlen=scalar @a;
    shift @a;
    $MaxQ=$Qlen if $Qlen>$MaxQ;
    my $Sum=0;
    $Sum += $_ for @a;
    if ($Sum>0) {
        $_ /= $Sum for @a;
    }
    $Dat{$Qin}=\@a;
}

for my $i (sort {$a<=>$b} keys %Dat) {
    for my $j (2..$MaxQ) {
        print "$i\t$j\t$Dat{$i}->[$j-2]\n"
    }
    print "\n";
}

