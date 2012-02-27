#!/bin/env perl
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
our $opts='r:o:l:p:s:c:b';
our($opt_o, $opt_r, $opt_l, $opt_p, $opt_s, $opt_c, $opt_b);

#our $desc='';
our $help=<<EOH;
\t-p type of input files {(auto),sam,soap,fq} [ONLY pIRS Generated fq]
\t-r ref fasta file (./ref/human.fa) [.{gz,bz2} is OK]
\t-s trim SNP positions from (<filename>) in format /^ChrID\\tPos/
\t-l read length of reads (100)
\t-o output prefix (./matrix).{count,ratio}.matrix and .{stat,info}
\t-c ChrID list (./chrtouse)
\t-b No pause for batch runs
For gzipped files, use zcat and pipe(|).
EOH
our $ARG_DESC='{sam,soap}pe_files';

ShowHelp();
$opt_r='./ref/human.fa' if ! $opt_r;
$opt_p='auto' if ! $opt_p;
$opt_o='./matrix' if ! $opt_o;
#$opt_c='./chrtouse' if ! $opt_c;
$opt_l=100 if ! $opt_l;
die "[x]-r $opt_r not exists !\n" unless -f $opt_r;
if ($opt_s) {die "[x]-s $opt_s not exists !\n" unless -f $opt_s;}

print STDERR "From [@ARGV]($opt_p) of [$opt_l] with [$opt_r] to [$opt_o]\n";
print STDERR "ChrID list:[$opt_c]\n" if $opt_c;
print STDERR "SNP skipping list:[$opt_s]\n" if $opt_s;
unless ($opt_b) {print STDERR "Wait 3 seconds to continue...\n"; sleep 3;}

#my $start_time = [gettimeofday];
#BEGIN

my %Genome;
if ($opt_c) {
    open C,'<',$opt_c or die "Error: $!\n";
    while(<C>){
        chomp;
        ++$Genome{$_};
    }
    close C;
}
warn "[!]Reading Reference Genome:\n";
if ($opt_r =~ /.bz2$/) {
    open( GENOME,"-|","bzip2 -dc $opt_r") or die "Error opening $opt_r: $!\n";
} elsif ($opt_r =~ /.gz$/) {
 	open( GENOME,"-|","gzip -dc $opt_r") or die "Error opening $opt_r: $!\n";
} else {open( GENOME,"<",$opt_r) or die "Error opening $opt_r: $!\n";}
#open GENOME,'<',$opt_r or die "Error: $!\n";
while (<GENOME>) {
    s/^>//;
	/^(\S+)/ or next;
	my $seqname = $1;
    print STDERR " >$seqname ...";
	$/=">";
	my $genome=<GENOME>;
	chomp $genome;
	$genome=~s/\s//g;
	$/="\n";
	if ((!$opt_c) or exists $Genome{$seqname}) {
        $Genome{$seqname}=$genome;
        print STDERR "\b\b\b",length $Genome{$seqname},".\n";
    } else {print STDERR "\b\b\b",length $genome,", skipped.\n";}
	$genome='';
}
close GENOME;
if ($opt_s) {
    print STDERR "[!]Reading SNP: ";
    open SNP,'<',$opt_s or die "Error: $!\n";
    while(<SNP>) {
        chomp;
        my ($chr,$pos)=split /\s+/;
        substr $Genome{$chr},$pos-1,1,'x' if exists $Genome{$chr};
    }
    close SNP;
    print STDERR "done.\n";
}
###
#print ">$_\n$Genome{$_}\n\n" for sort keys %Genome;
###
sub getBases($$$) {
    my ($chr,$start,$len)=@_;
    return substr $Genome{$chr},$start-1,$len;
}

my $READLEN=$opt_l;
my $MaxQ=40;
my $MinQ=3;
my $MisBase=0;
my $CountGridSampledOK_MinValue = 100;
my @SuggestGSPercent = (90, 5, 3);

=pod
http://en.wikipedia.org/wiki/FASTQ_format#Encoding
 S - Sanger        Phred+33,  raw reads typically (0, 40)
 X - Solexa        Solexa+64, raw reads typically (-5, 40)
 I - Illumina 1.3+ Phred+64,  raw reads typically (0, 40)
 J - Illumina 1.5+ Phred+64,  raw reads typically (3, 40)
    with 0=unused, 1=unused, 2=Read Segment Quality Control Indicator (bold) 
    (Note: See discussion above).
 L - Illumina 1.8+ Phred+33,  raw reads typically (0, 41)

Only J supported now.
=cut
my ($TotalBase,$TotalReads,%BaseCountTypeRef);
my ($mapBase,$mapReads,$QBbase,$QBmis)=(0,0,0,0);
my $Qascii=33;  # Sam 33, Soap 64.
my %Stat;   # $Stat{Ref}{Cycle}{Read}{Quality}
my %MarkovStat;   # $Stat{Ref}{Cycle}{pre-Q}{Read}{Quality}
my %QTrans;   # $QTrans{Cycle}{pre-Q}{Quality}
my %PlotReadsQavgHist;   # $PlotReadsQavgHist{Read1,2}{Q(round to 0.5)}=count, -1 => Sum (Illumina Q>0, but maybe 0 after round)
my %statQmkv;	# {$Len}->{PreQ_Avg}{Q}
my $QmkvMaxLen=3;

sub statRead($$$$$) {
    my ($ref,$isReverse,$read,$Qstr,$cyclestart)=@_;
    if ($isReverse) {
        $ref =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
        $read =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
		$ref = scalar reverse $ref;
		$read = scalar reverse $read;
		$Qstr = scalar reverse $Qstr;
    }
#	doTheStat($ref,$read,$Qstr,$cyclestart);
#}
#sub doTheStat($$$$) {
#    my ($ref,$read,$Qstr,$cyclestart)=@_;
    my $PEpos=-1;
    my $QBflag=0;
    my $lastQ=-1;
	my $SumQ=0;
	my @Qvalues;
    for (my $i=0;$i<$READLEN;$i++) {
        my $refBase=substr $ref,$i,1 or return;
        my $QstrSingle=substr $Qstr,$i,1;
        my $Qval=ord($QstrSingle)-$Qascii;
		push @Qvalues,$Qval;
		next unless $refBase =~ /^[ATCG]$/;
        my $readBase=substr $read,$i,1;
        next if $readBase eq 'N';
        if ($MaxQ<$Qval) {
            $MaxQ=$Qval;
            warn "[!] Qval=$Qval($QstrSingle) > 40 found. Remember to add -I at bwa aln for Illumina reads !\n";
        }
		$PEpos=$cyclestart+$i;
		if ($lastQ != -1) {	# 1st cycle skipped; 1st can be 'N', so not $PEpos > 1
			++$MarkovStat{$refBase}{$PEpos}{$lastQ}{$readBase}{$Qval};
			++$QTrans{$PEpos}{$lastQ}{$Qval};
		}
		$lastQ = $Qval;
		$SumQ += $Qval;
        ++$Stat{$refBase}{$PEpos}{$readBase}{$Qval};
        if ($Qval <= 2) {
            $QBflag = 1;
            ++$QBbase;
			$MinQ = $Qval if $MinQ > $Qval;	# well, if v1.8+, $MinQ can touch 0.
        }
        if ($refBase ne $readBase) {
            ++$MisBase;
            $QBmis += $QBflag;
        }
        ++$BaseCountTypeRef{$refBase};
        ++$TotalBase;
#print "$isReverse {$refBase}{$PEpos}{$readBase}{$Qval} ",($refBase eq $readBase)?'=':'x',"\n";
    }
	my $Read_num;
	if ($PEpos > $READLEN) {
		$Read_num = 2;
	} else {
		$Read_num = 1;
	}
    ++$TotalReads unless $PEpos==-1;
	++$PlotReadsQavgHist{$Read_num}{int(2*$SumQ/$READLEN)/2};	# 0 for [0,0.5]
	++$PlotReadsQavgHist{$Read_num}{-1};

	for my $qlen (1..$QmkvMaxLen) {
		my $pQavg=0;
		for (my $i=0;$i < $READLEN - $qlen + 1;$i++) {
			$pQavg += $Qvalues[$_] for ($i .. $i+$qlen-1);
			++$statQmkv{$qlen}{$Qvalues[$i]}->[int(2*$pQavg/$qlen)];
		}
	}
}

my $type;
if ($opt_p eq 'sam') {
    $type = 'sam';
} elsif ($opt_p eq 'soap') {
    $type = 'soap';
    $Qascii = 64;
} elsif ($opt_p eq 'fq') {
    $type = 'fq';
} else {
    chomp($_=<>) or die "[x]Empty input !\n";
    if (/^@/) {
        if (/^@\w+\t/) {
            $type = 'sam';
        } else {
            $type = 'fq';
            # @read_800_22/1 ChrID 2 129025 70 786 32,C;70,A;
        }
    } else {
        $type = 'soap';
        $Qascii = 64;
        my @read1=split /\t/;
        chomp($_=<>) or die "[x]Input too short !\n";
        my @read2=split /\t/;
        die "[x]Not PE soap file.\n" if $read1[0] ne $read2[0];
        $mapBase += $read1[5]+$read2[5];
        $mapReads +=2;
        goto LABEL unless exists $Genome{$read1[7]};
        goto LABEL unless $read1[3] == 1 and $read2[3] == 1;  # Hit==1
        goto LABEL if $read1[9] > 100 or $read2[9] > 100; # No Indel
        goto LABEL unless $read1[5] == $READLEN;
        goto LABEL unless $read2[5] == $READLEN;
        goto LABEL unless $read1[7] eq $read2[7];   # same Reference sequence NAME
        my $ref1=uc getBases($read1[7],$read1[8],$READLEN) or print join("\t",@read1),"\n";
        my $ref2=uc getBases($read2[7],$read2[8],$READLEN) or print join("\t",@read2),"\n";
        #my ($QNAME,$Seq,$Qual,$Hit,$a/b,$Len,$Strand,$Chr,$Pos,$Type,$SMID,$CIAGR,$Etc)=@read1;
        #       0     1    2     3    4    5     6      7    8    9     10     11   12
        statRead($ref1,$read1[6] eq '-',$read1[1],$read1[2],1);
        statRead($ref2,$read2[6] eq '-',$read2[1],$read2[2],1+$READLEN);
    }
}
LABEL:
print STDERR "[!]Input file type is [$type].\n";
my $start_time = [gettimeofday];

#my ($RL1,$RL2)=(0,0);
if ($type eq 'sam') {
    while (<>) {
        next if /^@\w\w\t\w\w:/;
        chomp;
        my @read1=split /\t/;
        chomp($_=<>) or last;
        my @read2=split /\t/;
    #print join("\t",@read1),"\n-",join("\t",@read2),"\n";
        die "[x]Not PE sam file.\n" if $read1[0] ne $read2[0];
        ++$mapReads if $read1[2] ne '*';
        ++$mapReads if $read2[2] ne '*';
        next unless exists $Genome{$read1[2]};
        next unless ($read1[1] & 3) == 3;  # paired + mapped in a proper pair
        next if $read1[1] >= 256;   # not primary || QC failure || optical or PCR duplicate
        next unless ($read2[1] & 3) == 3;
        next if $read2[1] >= 256;
        #$RL1=length($read1[9]);
        #$RL2=length($read2[9]);
        next unless $read1[5] =~ /^(\d+)M$/;
        next unless $1 == $READLEN;
        next unless $read2[5] =~ /^(\d+)M$/;
        next unless $1 == $READLEN;
        next unless $read1[6] eq '=';   # same Reference sequence NAME
        next unless $read2[6] eq '=';
        next if $read1[11] eq 'XT:A:R'; # Type: Unique/Repeat/N/Mate-sw, N not found.
        next if $read2[11] eq 'XT:A:R';
        my $ref1=uc getBases($read1[2],$read1[3],$READLEN) or print join("\t",@read1),"\n";
        my $ref2=uc getBases($read2[2],$read2[3],$READLEN) or print join("\t",@read2),"\n";
        #my ($QNAME,$FLAG,$RNAME,$POS,$MAPQ,$CIAGR,$MRNM,$MPOS,$ISIZE,$SEQ,$QUAL,$OPT)=@read1;
        #       0      1    2       3   4       5   6       7     8     9    10    11
        statRead($ref1,$read1[1] & 16,$read1[9],$read1[10],1);
        statRead($ref2,$read2[1] & 16,$read2[9],$read2[10],1+$READLEN);
    }
} else {
    while (<>) {
        #next if /^@\w\w\t\w\w:/;
        chomp;
        my @read1=split /\t/;
        chomp($_=<>) or last;
        my @read2=split /\t/;
    #print join("\t",@read1),"\n-",join("\t",@read2),"\n";
        die "[x]Not PE soap file.\n" if $read1[0] ne $read2[0];
        $mapBase += $read1[5]+$read2[5];
        $mapReads +=2;
        next unless exists $Genome{$read1[7]};
        next unless $read1[3] == 1 and $read2[3] == 1;  # Hit==1
        next if $read1[9] > 100 or $read2[9] > 100; # No Indel
        next unless $read1[5] == $READLEN;
        next unless $read2[5] == $READLEN;
        next unless $read1[7] eq $read2[7];   # same Reference sequence NAME
        my $ref1=uc getBases($read1[7],$read1[8],$READLEN) or print join("\t",@read1),"\n";
        my $ref2=uc getBases($read2[7],$read2[8],$READLEN) or print join("\t",@read2),"\n";
        #my ($QNAME,$Seq,$Qual,$Hit,$a/b,$Len,$Strand,$Chr,$Pos,$Type,$SMID,$CIAGR,$Etc)=@read1;
        #       0     1    2     3    4    5     6      7    8    9     10     11   12
        statRead($ref1,$read1[6] eq '-',$read1[1],$read1[2],1);
        statRead($ref2,$read2[6] eq '-',$read2[1],$read2[2],1+$READLEN);
    }
}

open OA,'>',$opt_o.'.count.matrix' or die "Error: $!\n";
open OB,'>',$opt_o.'.ratio.matrix' or die "Error: $!\n";
open OC,'>',$opt_o.'.stat' or die "Error: $!\n";
open OI,'>',$opt_o.'.info' or die "Error: $!\n";
my (%CountGridOK,%CountGridPoor,%CountGridZero,%CountGridSampled);

sub toCountGridSampled($$) {
	my ($c,$type)=@_;
	if ($c==0) {
		++$CountGridZero{$type};
	} elsif ($c < $CountGridSampledOK_MinValue) {
		++$CountGridPoor{$type};
	} else {
		++$CountGridOK{$type};
	}
	++$CountGridSampled{$type}{$c};
}

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
my $Qcount=$MaxQ+1;
$TotalBase=-1 unless $TotalBase;
my $MisRate=100*$MisBase/$TotalBase;
print OI "Profile Name: ${opt_o}.(count|ratio).matrix, Generate @ $date by ${user}$mail .\n\n";
$tmp="#Generate @ $date by ${user}$mail
#Input [$type] file of mapped Reads: $mapReads , mapped Bases $mapBase (no base stat for sam files)
#Total statistical Bases: $TotalBase , Reads: $TotalReads of ReadLength $READLEN
#Dimensions: Ref_base_number 4, Cycle_number $Cycle, Seq_base_number 4, Quality_number $Qcount
#Mismatch_base: $MisBase, Mismatch_rate: $MisRate %
#QB_Bases: $QBbase, QB_Mismatches: $QBmis (bases with quality <= 2)
#Reference Base Ratio in reads: ";
my @BaseOrder=sort qw{A T C G}; # keys %BaseCountTypeRef;
my $RBR;
for (@BaseOrder) {
    $RBR .= $_.' '. int(0.5+100*1000*$BaseCountTypeRef{$_}/$TotalBase)/1000 .' %;   ';
}
$tmp .= $RBR;
my @BaseQ;
for my $base (@BaseOrder) {
    push @BaseQ,"$base-$_" for (0..$MaxQ);
}
$tmp .= "\n
[Info]
Date = $date
User = ${user}$mail
FileType = $type
ReadLength = $READLEN
Ref_base_number = 4
Cycle_number = $Cycle
Seq_base_number = 4
Quality_number = $Qcount
Mismatch_rate = $MisRate %
Reference_Base_Ratio = $RBR
<<END

[Stat]
MappedReads = $mapReads
MappedBases = $mapBase
UsedReads = $TotalReads
UsedBase = $TotalBase
MismatchBase = $MisBase
QB_Bases = $QBbase
QB_Mismatches = $QBmis
<<END

[DistMatrix]
#".join("\t",'Ref','Cycle',@BaseQ);
print OA $tmp;
print OB $tmp;
print OA "\tRowSum\n";
print OB "\n";
my ($count,$countsum);
for my $ref (@BaseOrder) {
    #print OA "\n";
    for my $cycle (1..(2*$READLEN)) {
        $tmp="$ref\t$cycle\t";
        print OA $tmp; print OB $tmp;
        my (@Counts,@Rates)=();
        for my $base (@BaseOrder) {
            for my $q (0..$MaxQ) {
                if (exists $Stat{$ref}{$cycle} and exists $Stat{$ref}{$cycle}{$base} and exists $Stat{$ref}{$cycle}{$base}{$q}) {
                    $count=$Stat{$ref}{$cycle}{$base}{$q};
                } else {$count=0;}
                push @Counts,$count;
                &toCountGridSampled($count,'4D') if $q >= $MinQ;
            }
        }
        $countsum=0;
        $countsum += $_ for @Counts;
        $countsum=-1 if $countsum==0;
        push @Rates,$_/$countsum for @Counts;
        print OA join("\t",@Counts,$countsum),"\n";
        print OB join("\t",@Rates),"\n";
    }
}
print OA "<<END\n\n";

print OA "[QTransMatrix]\n#",join("\t",'Cycle','pre1_Q',0..$MaxQ),"\tRowSum\n";
#for my $ref (@BaseOrder) {
    for my $cycle (2..(2*$READLEN)) {
		next if $cycle == 1 + $READLEN;	# the first cycle is always 0
		for my $preQ ($MinQ..$MaxQ) {
			print OA "$cycle\t$preQ\t";
			my @Counts=();
			for my $q (0..$MaxQ) {
				if (exists $QTrans{$cycle} and exists $QTrans{$cycle}{$preQ} and exists $QTrans{$cycle}{$preQ}{$q}) {
					$count=$QTrans{$cycle}{$preQ}{$q};
				} else {$count=0;}
				push @Counts,$count;
				&toCountGridSampled($count,'QT') if $q >= $MinQ;
			}
			$countsum=0;
			$countsum += $_ for @Counts;
			print OA join("\t",@Counts,$countsum),"\n";
		}
    }
#}
print OA "<<END\n";
print OB "<<END\n";
close OA;
close OB;

print OC "[5Dmatrix]\nType = 3d\n#",join("\t",'Ref','Cycle','pre1_Q',@BaseQ),"\tRowSum\n";
for my $ref (@BaseOrder) {
    for my $cycle (1..(2*$READLEN)) {
		for my $preQ ($MinQ..$MaxQ) {
			print OC "$ref\t$cycle\t$preQ\t";
			my @Counts=();
			for my $base (@BaseOrder) {
				for my $q (0..$MaxQ) {
					if (exists $MarkovStat{$ref}{$cycle} and exists $MarkovStat{$ref}{$cycle}{$preQ} and exists $MarkovStat{$ref}{$cycle}{$preQ}{$base} and exists $MarkovStat{$ref}{$cycle}{$preQ}{$base}{$q}) {
						$count=$MarkovStat{$ref}{$cycle}{$preQ}{$base}{$q};
					} else {$count=0;}
					push @Counts,$count;
					&toCountGridSampled($count,'5D') if $q >= $MinQ;
				}
			}
			$countsum=0;
			$countsum += $_ for @Counts;
			print OC join("\t",@Counts,$countsum),"\n";
		}
    }
}
print OC "<<END\n";

print OC "\n[AvgQonReads]
#Total Quality values: $PlotReadsQavgHist{1}{-1}, $PlotReadsQavgHist{2}{-1}
#Q\tRead_1\tRead_2\tRatio_1\tRatio_2\n";
for my $q (0..2*$MaxQ) {
	$PlotReadsQavgHist{1}{$q/2} = 0 unless exists $PlotReadsQavgHist{1}{$q/2};
	$PlotReadsQavgHist{2}{$q/2} = 0 unless exists $PlotReadsQavgHist{2}{$q/2};
	print OC join("\t",$q/2,$PlotReadsQavgHist{1}{$q/2},$PlotReadsQavgHist{2}{$q/2},
		$PlotReadsQavgHist{1}{$q/2}/$PlotReadsQavgHist{1}{-1},
		$PlotReadsQavgHist{2}{$q/2}/$PlotReadsQavgHist{2}{-1}),"\n";
}
print OC "<<END\n";

print OC "\n[QtransStat]\n";
my @t;
push @t,$_/2 for (2*2..2*$MaxQ);
print OC join("\t",'Q','preLen','preQmean',@t),"\n";
for my $qlen (sort {$a<=>$b} keys %statQmkv) {
	for my $currQ (sort {$a<=>$b} keys %{$statQmkv{$qlen}}) {
		print OC "$currQ\t$qlen\t";
		@t=();
		my ($sumT,$sumQ)=(0,0);
		for (2*2..2*$MaxQ) {
			if (defined $statQmkv{$qlen}{$currQ}->[$_/2]) {
				push @t,$statQmkv{$qlen}{$currQ}->[$_/2];
				$sumT += $statQmkv{$qlen}{$currQ}->[$_/2];
				$sumQ += $currQ * ($statQmkv{$qlen}{$currQ}->[$_/2]);
			} else {
				push @t,'-';
			}
		}
		#$sumT /= 2*$MaxQ - 2*$MinQ +1;	# This is just average count, >_<
		if ($sumT != 0) {
			print OC join("\t",int(0.5+10000*$sumQ/$sumT)/10000,@t),"\n";;
		} else {
			print OC join("\t",'-',@t),"\n";
		}
	}
}
print OC "<<END\n";

close OC;

sub getValueNoNULL($) {
	return $_[0] if defined $_[0];
	return 0;
}

my $TotalGrid = 1;
my ($CountGridOK,$CountGridPoor,$CountGridZero);
for $type (sort keys %CountGridSampled) {
	$CountGridOK = getValueNoNULL($CountGridOK{$type});
	$CountGridPoor = getValueNoNULL($CountGridPoor{$type});
	$CountGridZero = getValueNoNULL($CountGridZero{$type});
	#($CountGridOK,$CountGridPoor,$CountGridZero)=($CountGridOK{$type},$CountGridPoor{$type},$CountGridZero{$type});
	if ($type eq '4D') {
		$TotalGrid = 16 * 2 * $READLEN * ($MaxQ -$MinQ + 1);
	} elsif ($type eq '5D') {
		$TotalGrid = $CountGridOK+$CountGridPoor+$CountGridZero;
	}
	print OI "[$type]\nTotal_Grid: $TotalGrid
	Sampling_Thresholds: $CountGridSampledOK_MinValue
	Sampled_Enough: $CountGridOK (",int(10000*$CountGridOK/$TotalGrid)/100,"%), should be > $SuggestGSPercent[0] %
	Sampled_Poor: $CountGridPoor (",int(10000*$CountGridPoor/$TotalGrid)/100,"%) should be < $SuggestGSPercent[1] %
	Empty_Grid: $CountGridZero (",int(10000*$CountGridZero/$TotalGrid)/100,"%) should be < $SuggestGSPercent[2] %
	\n";
	if ($CountGridPoor) {
		print OI "Poor Sampling Histogram:\n";
		for (1 .. $CountGridSampledOK_MinValue) {
			next unless exists $CountGridSampled{$type}{$_};
			print OI "$_\t$CountGridSampled{$type}{$_}\n";
		}
	}
	print OI "\nYou may need to supply more data to get a better profile.\n" if 100 * $CountGridOK < $SuggestGSPercent[0];
	print OI "\n";
}
close OI;
#END
my $stop_time = [gettimeofday];

print STDERR "\nTime Elapsed:\t",tv_interval( $start_time, $stop_time )," second(s).\n";
__END__
zcat bwamask/mask110621_I263_FCB066DABXX_L8_HUMjrmRACDKAAPEI-3.sam.gz|head -n200|./matrixsam.pl -b 2>&1 |tee logerr.txt

3289m for Hg18 reference.

0.192234941 ms per line.
Thus, for a LANE of 216009064 lines, which is (108004507 sequences)x2+50, 11.5345804648626 hours needed.

./baseCalling_Matrix_calculator.pl -c chrtouse -b -o test2 t.sam
