#!/usr/bin/env perl
=pod
Author: Hu Xuesong @ BGI <huxuesong@genomics.org.cn>
Version: 0.1.1 @ 20110819
=cut
use strict;
use warnings;
use Time::HiRes qw ( gettimeofday tv_interval );
use Getopt::Std;

$Getopt::Std::STANDARD_HELP_VERSION=1;
sub main::ShowHelp() {
	if (@main::ARGV == 0) {
		&main::VERSION_MESSAGE();
		&main::HELP_MESSAGE();
		die "\n";
	}
	getopts($main::opts);
}

sub main::HELP_MESSAGE() {
	$main::help =~ s|\[|[\033[0;0m|g;
	$main::help =~ s|\]|\033[32;1m]|g;
	$main::help =~ s|\(|(\033[0;1m|g;
	$main::help =~ s|\)|\033[32;1m)|g;
	$main::help =~ s|:(\s*\n?\s*)(\S)|:$1\033[0;1m$2|g;
	$main::help =~ s|\\\[\033\[0;0m|[|g;
	$main::help =~ s|\\\033\[32;1m\]|]|g;
	$main::help =~ s|\\\(\033\[0;1m|(|g;
	$main::help =~ s|\\\033\[32;1m\)|)|g;
	$main::help =~ s|\\:(\s*\n?\s*)\033\[0;1m|:$1|g;
	$main::help =~ s|\n|\033[32;1m\n|g;
	$main::ARG_DESC='[PROGRAM_ARG1 ...]' unless $main::ARG_DESC;
	print STDERR <<EOH;
\nUsage: \033[0;1m$0\033[0;0m [-OPTIONS [-MORE_OPTIONS]] [--] $main::ARG_DESC

The following single-character options are accepted:
\033[32;1m$main::help\033[0;0mOptions may be merged together.  -- stops processing of options.
Space is not required between options and their arguments.
EOH
}
sub main::VERSION_MESSAGE() {
	my $perlv = $];
	$perlv = sprintf "%vd", $^V if $] >= 5.006;
	my $ver = sprintf "%vd", $main::VERSION;
	my ($scr) = ($0 =~ m,([^/\\]+)$,);
	if ($main::desc) {
		print STDERR <<EOH;
\033[32;1m$main::desc\033[0;0m ($scr) version \033[0;1m$ver\033[0;0m,
 running under Perl version $perlv.
EOH
	} else {
		print STDERR <<EOH;
\033[32;1m$scr\033[0;0m version \033[0;1m$ver\033[0;0m, running under Perl version $perlv.
EOH
	}
}

$main::VERSION=0.1.1;
our $opts='o:b';
our($opt_o, $opt_b);

#our $desc='';
our $help=<<EOH;
\t-o output prefix (./allmatrix).{count,ratio}.matrix
\t-b No pause for batch runs
EOH
our $ARG_DESC='matrix_count_files';

ShowHelp();
$opt_o='./allmatrix' if ! $opt_o;

print STDERR "From [@ARGV] to [$opt_o]\n";
unless ($opt_b) {print STDERR 'press [Enter] to continue...'; <STDIN>;}

my $start_time = [gettimeofday];
#BEGIN
my $READLEN=0;
my $Qcount=41;
my ($TotalReads,$TotalBase,$MisBase,%BaseCountTypeRef)=(0,0,0);
my ($mapBase,$mapReads,$QBbase,$QBmis)=(0,0,0,0);
my $type='N/A';
my %Stat;   # $Stat{Ref}{Cycle}{Read-Quality}
my @BQHeader;
while (<>) {
    if (/^#Input \[(\w+)\] file of mapped Reads: (\d+) , mapped Bases (\d+)/) {
        $type=$1 if $type ne 'sam';
        $mapReads += $2;
        $mapBase += $3;# if $type ne 'sam';
    }
    if (/^#Total statistical Bases: (\d+) , Reads: (\d+) of ReadLength (\d+)$/) {
        print " >$ARGV|   Reads:$2   ReadLen:$3   Bases:$1\n";
        $TotalReads += $2;
        $READLEN = $3 if $READLEN < $3;
    }
    if (/^#Ref\tCycle\t/) {
        #s/^#//;
        chomp;
        (undef,undef,@BQHeader)=split /\t/;
        pop @BQHeader if $BQHeader[-1] eq 'RowSum';
    }
    if (/^#Dimensions:.+?Quality_number (\d+)$/) {
        $Qcount = $1 if $Qcount<$1;
    }
    if (/^#Mismatch_base: (\d+)/) {
        $MisBase += $1;
    }
    if (/^#QB_Bases: (\d+), QB_Mismatches: (\d+)/) {
        $QBbase += $1;
        $QBmis += $2;
    }
    next if /^#/;
    next if /^$/;
    chomp;
    my ($ref,$cycle,@BQ)=split /\t/;
    #print "$ref,$cycle,@BQ\n";
    next unless $ref =~ /^[ATCG]$/;
    #die "[$_]\n$ref,$cycle,[@BQ]\n$#BQ < $#BQHeader " if $#BQ < $#BQHeader;
    for my $key (@BQHeader) {
        my $value=shift @BQ;
        $Stat{$ref}{$cycle}{$key}+=$value;
        $BaseCountTypeRef{$ref}+=$value;
        $TotalBase+=$value;
        #print "{$ref}{$cycle}{$key}$value\n";
    }
}
#print $TotalReads,"\t",$READLEN,"\n";
#print join("\t",@BQHeader),"\n";
open OA,'>',$opt_o.'.count.matrix' or die "Error: $!\n";
open OB,'>',$opt_o.'.ratio.matrix' or die "Error: $!\n";
my $tmp;
chomp(my $user=`id -nru`);
@ARGV=('/etc/passwd');
chomp(my @passwd=<>);
($_)=grep /$user/,@passwd;
$tmp=(split /:/)[4];
my $mail='';
$mail=" <$tmp>" unless $tmp =~ /^\s*$/;
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
my $date=sprintf "%02d:%02d:%02d,%4d-%02d-%02d",$hour,$min,$sec,$year+1900,$mon+1,$mday;
my $Cycle=2*$READLEN;
my $MisRate=100*$MisBase/$TotalBase;
$tmp="#Generate @ $date by ${user}$mail
#Input [$type] file of Reads: $mapReads , Bases $mapBase (no base stat for sam files)
#Total statistical Bases: $TotalBase , Reads: $TotalReads of ReadLength $READLEN
#Dimensions: Ref_base_number 4, Cycle_number $Cycle, Seq_base_number 4, Quality_number $Qcount
#Mismatch_base: $MisBase, Mismatch_rate: $MisRate %
#QB_Bases: $QBbase, QB_Mismatches: $QBmis (bases with quality <= 2)
#Reference Base Ratio in reads: ";
my @BaseOrder=sort keys %BaseCountTypeRef;  # qw{A T C G};
for (@BaseOrder) {
    $tmp .= $_.' '. int(0.5+100*1000*$BaseCountTypeRef{$_}/$TotalBase)/1000 .' %;   ';
}

$tmp .= "\n#".join("\t",'Ref','Cycle',@BQHeader);
print OA $tmp;
print OB $tmp;
print OA "\tRowSum";
print OB "\n";
my ($count,$countsum);
for my $ref (@BaseOrder) {
    print OA "\n";
    for my $cycle (sort {$a<=>$b} keys %{$Stat{$ref}}) {
        $tmp="$ref\t$cycle\t";
        print OA $tmp; print OB $tmp;
        my (@Counts,@Rates)=();
        for my $bq (@BQHeader) {
            push @Counts,$Stat{$ref}{$cycle}{$bq};
        }
        #print "[",join("|",@Counts),"\n";
        $countsum=0;
        $countsum += $_ for @Counts;
        push @Rates,$_/$countsum for @Counts;
        print OA join("\t",@Counts,$countsum),"\n";
        print OB join("\t",@Rates),"\n";
    }
}
close OA;
close OB;

#END
my $stop_time = [gettimeofday];

print STDERR "\nTime Elapsed:\t",tv_interval( $start_time, $stop_time )," second(s).\n";
__END__

