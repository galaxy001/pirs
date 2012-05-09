#!/usr/bin/env perl
=pod
Author: Hu Xuesong @ BGI <huxuesong@genomics.org.cn>
Version: 1.0.0 @ 20110803
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

$main::VERSION=1.0.0;
our $opts='o:b';
our($opt_o, $opt_b);

#our $desc='';
our $help=<<EOH;
\t-o output prefix (./stat).{info,insD}
\t-b No pause for batch runs
For SOAP, both soap/single files and STDERR dump are needed.
For BWA , only insert size stat. is done currently.
EOH
our $ARG_DESC='soap/sam files{,.gz,.bz2}';

ShowHelp();
$opt_o='./stat' if ! $opt_o;
die "[x]No input files found !\n" unless @ARGV;
#die "[!]Max 252 files supported.\n" if @ARGV>252;

print STDERR "From [@ARGV] to [$opt_o]\n";
unless ($opt_b) {print STDERR 'press [Enter] to continue...'; <STDIN>;}

my $start_time = [gettimeofday];
#BEGIN
sub combineC($) {
	my $href=$_[0];
	if ($href and %$href) {
		my (@str,$m);
		$m = (sort {$a<=>$b} keys %$href)[-1];
		for (1..$m) {
			#$$href{$_} += 0;
			push @str,join(':',$_,$$href{$_}||0);
		}
		return \join(',',@str);
	} else {return \'.';}
}

sub combineJ($) {
	my $href=$_[0];
	if ($href and %$href) {
		my @str;
		for (sort {$a<=>$b} keys %$href) {
			push @str,join(':',$_,$$href{$_});
		}
		return \join(',',@str);
	} else {return \'.';}
}

sub getRealpos($$$$) {
	my ($len,$strand,$realpos,$trim)=@_;
	if ($strand eq '-') {	# Negative
		$realpos += $len;	# should be $len-1. So, starting 0. (+ & -, never meets.)
		if ($trim =~ /(\d+)S$/) {
			$realpos += $1;	# $1 only reset after next /()/
		}
	} elsif ($strand eq '+') {	# Positive
		if ($trim =~ /^(\d+)S/) {
			$realpos -= $1;
		}
	} else {
		$realpos=-1;
	}
	return $realpos;
}

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
sub checkfiletype($) {
    my $fhc=$_[0];
    my $type;
    my $head=<$fhc>;
    if ($head =~ /^\[bwa_sai2sam_pe_core\]/) {
        $type='samlog';
    } elsif ($head =~ /^\@SQ/) {
        $type='sam';
    } elsif ($head =~ /[ATCG]/) {
        $type='soap';
    } elsif ($head =~ /^$/) {
=pod
00000000  0a 42 65 67 69 6e 20 50  72 6f 67 72 61 6d 20 53  |.Begin Program S|
00000010  4f 41 50 61 6c 69 67 6e  65 72 2f 73 6f 61 70 32  |OAPaligner/soap2|
00000020  0a 46 72 69 20 4a 75 6c  20 31 35 20 30 39 3a 32  |.Fri Jul 15 09:2|
00000030  39 3a 34 32 20 32 30 31  31 0a 52 65 66 65 72 65  |9:42 2011.Refere|
00000040  6e 63 65 3a 20 2e 2f 72  65 66 2f 68 75 6d 61 6e  |nce: ./ref/human|
00000050  2e 66 61 2e 69 6e 64 65  78 0a 51 75 65 72 79 20  |.fa.index.Query |
00000060  46 69 6c 65 20 61 3a 20  2e 2f 30 6e 6f 6d 61 73  |File a: ./0nomas|
=cut
        $type='soaplog';
    }
    return $type;
}

sub sumsoapdata($$) {
    my ($hit,$len,$chr,$types,$trim,$mistr)=@{$_[0]};
    my $dathref=$_[1];
    my ($BPOut,$ReadsOut,$MisSum,$TrimedBP,$TrimedReads,%Hit9r,%Hit9bp,%misMatch,%Indel)=(0,0,0,0,0);
    #my (%chrBPOut,%chrReadsOut,%chrMisSum,%chrTrimedBP,%chrTrimedReads,%chrHit9r,%chrHit9bp,%chrmisMatch,%chrIndel);
    my (@trims,$missed);
		$BPOut += $len;#	$chrBPOut{$chr} += $len;
		++$ReadsOut;#	++$chrReadsOut{$chr};
		$missed=$mistr=~tr/ATCGatcg//;
		++$misMatch{$missed};
		#++$chrmisMatch{$chr}{$missed};
		$MisSum += $missed;
		#$chrMisSum{$chr} += $missed;

		$hit=4 if $hit>4;	# max to count 3, then >=4. Ancient Chinese wisdom, and a bit more ...
		++$Hit9r{$hit};
		#++$chrHit9r{$chr}{$hit};
		$Hit9bp{$hit} += $len;
		#$chrHit9bp{$chr}{$hit} += $len;
=pod
		if ($types < 100) {
			#++$misMatch{$types};
			#++$chrmisMatch{$chr}{$types};
		} elsif ($types < 200) {	# '3S33M9D39M', '32M1D14M29S' exists ?
			++$Indel{$types-100};
			++$chrIndel{$chr}{$types-100};
		} else {
			++$Indel{200-$types};
			++$chrIndel{$chr}{200-$types};
		}
=cut
		if ($types > 200) {
			++$Indel{200-$types};
			#++$chrIndel{$chr}{200-$types};
		} elsif ($types > 100) {
			++$Indel{$types-100};
			#++$chrIndel{$chr}{$types-100};
		}

		@trims = $trim =~ /(\d+)S/;
		if (@trims) {
			++$TrimedReads;
			#++$chrTrimedReads{$chr};
			for ( @trims ) {
				$TrimedBP += $_;
				#$chrTrimedBP{$chr} += $_;
			}
		}
    $dathref->{'BPOut'} += $BPOut;
    $dathref->{'ReadsOut'} += $ReadsOut;
    $dathref->{'MisSum'} += $MisSum;
    $dathref->{'TrimedReads'} += $TrimedReads;
    $dathref->{'TrimedBP'} += $TrimedBP;
    $dathref->{'misMatch'}{$_} += $misMatch{$_} for keys %misMatch;
    $dathref->{'Indel'}{$_} += $Indel{$_} for keys %Indel;
    $dathref->{'Hit9r'}{$_} += $Hit9r{$_} for keys %Hit9r;
    $dathref->{'Hit9bp'}{$_} += $Hit9bp{$_} for keys %Hit9bp;
}
sub statsoap($) {
    my $fh=$_[0];
    my (%datsum);
    $datsum{'PEuniqPairs'}=0;
    my ($BadLines,$PESE,$pairs,$lastpos,$line1,$line2,$pp,$pn,$calins,%insD)=(0,'PE',0);
	while ($line1=<$fh>) {
#print "[$line1]\n";
		my ($id1, $n1, $len1, $f1, $chr1, $x1, $types1, $m1, $mistr1)
		 = (split "\t", $line1)[0,3,5,6,7,8,9,-2,-1];
		unless (defined $types1) {    # soap2 output always more than 10 columes.
		    ++$BadLines;
		    last;
		}
		&sumsoapdata([$n1, $len1, $chr1, $types1, $m1, $mistr1],\%datsum);
		$line2=<$fh>;
		last unless $line2;
		my ($id2, $n2, $len2, $f2, $chr2, $x2, $types2, $m2, $mistr2)
		 = (split "\t", $line2)[0,3,5,6,7,8,9,-2,-1];
		unless (defined $types2) {    # soap2 output always more than 10 columes.
		    ++$BadLines;
		    last;
		}
			#($soapid,$hit,$len,$strand,$chr,$pos,$trim) = (split(/\t/))[0,3,5,6,7,8,-2];
			#        ($hit,$len,        $chr,$types,$trim,$mistr) = @lines[3,5,7,9,-2,-1];
		$id1 =~ s/\/[12]$//;
		$id2 =~ s/\/[12]$//;
		&sumsoapdata([$n2, $len2, $chr2, $types2, $m2, $mistr2],\%datsum);
		if (($PESE eq 'SE') or ($id1 ne $id2)){	# single
			$PESE='SE';
			next;
		}
		next if $n1+$n2>2 or $chr1 ne $chr2;
=pod
		if ($f1 eq '+') {
			($pp,$pn)=($x1,$x2);
		} else {
			($pp,$pn)=($x2,$x1);
		}
		if ($ins > 1500) {	# FR => +.pos < -.pos; RF => -.pos < +.pos
			next if $pp < $pn;
		} else { next if $pp > $pn; }
=cut
		++$pairs;
		$line1=&getRealpos($len1, $f1, $x1, $m1);	# $len,$strand,$realpos,$trim
		$line2=&getRealpos($len2, $f2, $x2, $m2);	# Well, $line{1,2} is recycled.
		$calins=abs($line1-$line2);	# -.starting=0
		++$insD{$calins};
	}
	$datsum{'BadLines'}=$BadLines;
	if ($PESE eq 'PE') {
    	$datsum{'PEuniqPairs'}=$pairs;
	}
	return [$PESE,[0,0,0],\%datsum,\%insD];
}
sub statsoaplog($) {
    my $fh=$_[0];
    my ($Pairs,$Paired,$Singled,$Reads,$Alignment)=(0,0,0,0,0);
    while(<$fh>) {
	    $Pairs = (split)[-2] if /^Total Pairs:/;
	    $Paired = (split)[1] if /^Paired:/;
	    $Singled = (split)[1] if /^Singled:/;
	    $Reads = (split)[-1] if /^Total Reads/;
	    $Alignment = (split)[1] if /^Alignment:/;
    }
    my @RET;
    if ($Reads) {
        @RET=('SE',[$Reads,0,$Alignment]);
    } else {
        @RET=('PE',[$Pairs*2,$Paired*2,$Singled]);
    }
    return \@RET;
=pod
Total Pairs: 34776407 PE
Paired:      17719335 (50.95%) PE
Singled:     30467634 (43.81%) SE

IN:	34776407 reads x 2 fq file = 69552814 reads from _1 & _2.fq
-o:	17719335 pairs = 17719335 x 2 =35438670 lines in .soap
-2:	30467634 lines in .single

17719335/34776407 = 0.50952172833726037310294878939046
30467634/69552814 = 0.43805034257851882168275750856033


Total Reads: 25
Alignment:   22 (88.00%)
=cut
}

sub sumsamdata($$) {
        #my ($QNAME,$FLAG,$RNAME,$POS,$MAPQ,$CIAGR,$MRNM,$MPOS,$ISIZE,$SEQ,$QUAL,$OPT)=@read1;
        #       0      1    2       3   4       5   6       7     8     9    10    11
    my ($FLAG,$RNAME,$POS,$CIAGR,$MRNM,$MPOS,$ISIZE)=@{$_[0]}[1,2,3,5,6,7,8];
    my $dathref=$_[1];
    my ($BPOut,$ReadsOut,$MisSum,$TrimedBP,$TrimedReads,%Hit9r,%Hit9bp,%misMatch,%Indel)=(0,0,0,0,0);
    my (@trims,@matches,@inserts,@deletes,$missed);
    my ($isPaired,$isSingled)=(0,0);
    unless ($CIAGR eq '*') {
        ++$ReadsOut;
    } else { return [0,0]; }
    if ($FLAG & 2) {
        $isPaired=1;
    } else {   #($FLAG & 8) and ($FLAG & 4 == 0)
        $isSingled=1;
    }
	@matches = $CIAGR =~ /(\d+)M/;
	for ( @matches ) {
		$BPOut += $_;
	}
	my ($Insertion,$Deletion)=(0,0);
	@inserts = $CIAGR =~ /(\d+)I/;
	for ( @inserts ) {
		$BPOut += $_;
		$Insertion += $_;
	}
	@deletes = $CIAGR =~ /(\d+)D/;
	for ( @deletes ) {
		$Deletion -= $_;
	}
	if ($Insertion and $Deletion) {
	    $Indel{$Insertion} += 0.5;
	    $Indel{$Deletion} += 0.5;
	} elsif ($Insertion) {
	    ++$Indel{$Insertion};
	} else {
	    ++$Indel{$Deletion};
	}
	my $Alternativehits='';
	for ( @{$_[0]} ) {
    	$Alternativehits = $1 if /^XA:Z:([\w,+-;]+)$/; #XA:Z:chrX,+1144912,100M,0;
	    next unless /^XM:i:(\d+)$/;
	    $MisSum += $1;
	}
	my @excigar=split(';',$Alternativehits);
	my $hit=1+scalar @excigar;
	$hit=4 if $hit>4;	# max to count 3, then >=4. Ancient Chinese wisdom, and a bit more ...
	++$Hit9r{$hit};
	$Hit9bp{$hit} += $BPOut;
	++$misMatch{$MisSum};
	@trims = $CIAGR =~ /(\d+)S/;
	if (@trims) {
		++$TrimedReads;
		for ( @trims ) {
			$TrimedBP += $_;
		}
	}
    $dathref->{'BPOut'} += $BPOut;
    $dathref->{'ReadsOut'} += $ReadsOut;
    $dathref->{'MisSum'} += $MisSum;
    $dathref->{'TrimedReads'} += $TrimedReads;
    $dathref->{'TrimedBP'} += $TrimedBP;
    $dathref->{'misMatch'}{$_} += $misMatch{$_} for keys %misMatch;
    $dathref->{'Indel'}{$_} += $Indel{$_} for keys %Indel;
    $dathref->{'Hit9r'}{$_} += $Hit9r{$_} for keys %Hit9r;
    $dathref->{'Hit9bp'}{$_} += $Hit9bp{$_} for keys %Hit9bp;
    return [$isPaired,$isSingled];
}
sub statsam($) {
    my $fh=$_[0];
    my ($PESE,$Reads,$toPaired,$toSingled,$BadLines,%datsum,%insD)=('PE',0,0,0,0);
    my ($Paired,$Singled)=(0,0);
    my ($line1,$line2,$calins,$sumret);
    while ($line1=<$fh>) {
        next if $line1=~/^@\w\w\t\w\w:/;
		my @read1=split /\t/, $line1;
		#print "[",join('|',@read1),"]\n";
		#print scalar @read1,"|\n";
		if (scalar @read1 < 11) {
		    ++$BadLines;
		    last;
		}
		++$Reads;
		$sumret=&sumsamdata(\@read1,\%datsum);
		$Paired += $sumret->[0];
		$Singled += $sumret->[1];
		$line2=<$fh>;
		last unless $line2;
		my @read2=split /\t/, $line2;
		if (scalar @read2 < 11) {
		    ++$BadLines;
		    last;
		}
		++$Reads;
		$sumret=&sumsamdata(\@read2,\%datsum);
		$Paired += $sumret->[0];
		$Singled += $sumret->[1];
		if (($PESE eq 'SE') or ($read1[0] ne $read2[0])){	# single
			$PESE='SE';
			next;
		}
		$calins=abs($read2[8]);
		++$insD{$calins} if $sumret->[0] and $calins < 1500;
    }
    $datsum{'BadLines'}=$BadLines;
    return [$PESE,[$Reads,$Paired,$Singled],\%datsum,\%insD];
}
my %dostat=(
    'sam'     => \&statsam,
    'soap'    => \&statsoap,
#    'samlog' => sub {},
    'soaplog' => \&statsoaplog,
);

sub sumintohash($$) {
    my ($inhref,$intohref)=@_;
    for my $key (keys %{$inhref}) {
        if (ref($$inhref{$key}) eq 'HASH') {
            $$intohref{$key}={} unless exists $$intohref{$key};
            &sumintohash($$inhref{$key},$$intohref{$key});
        } else {
            $$intohref{$key} += $$inhref{$key};
        }
    }
}

my $files=0;
my ($withPE,$InReads,$mapPair,$mapSingle,%DatSum,%InsD)=(0,0,0,0);
while($_=shift @ARGV) {
    next unless -f $_;
    ++$files;
    my $infile;
    $infile=openfile($_);
    my $type=checkfiletype($infile);
    close $infile;
    print STDERR "$files\t[$type] $_ ...";
    $infile=openfile($_);
    if (exists $dostat{$type}) {
        my $ret=&{$dostat{$type}}($infile); # [$PESE,[$Reads,$Paired*2,$Singled],\%datsum,\%insD]
        $InReads += $ret->[1]->[0];
        $mapPair += $ret->[1]->[1];
        $mapSingle += $ret->[1]->[2];
        if ($ret->[0] eq 'PE') {
            $withPE=1;
            &sumintohash($ret->[3],\%InsD) if $ret->[3];
        }
        &sumintohash($ret->[2],\%DatSum) if $ret->[2];
        #print "\n[";dump($ret);print "]\n";
    } else {
        print STDERR "\b\b\bskipped."
    }
    close $infile;
    print STDERR "\n";
}
warn "[!]Files Read: $files\n";

my ($n,$v,$avg,$std,$Lsd,$Rsd,$max_y,$max_x,$sum,$sum2)=(0);
if ($withPE) {
    open O,'>',"$opt_o.insD" or die "[x]Error opening $opt_o.insD: $!\n";
    for my $k (keys %InsD) {
	    $v=$InsD{$k};
	    $sum += $k * $v;
	    $n += $v;
	    $sum2 += $k*$k * $v;
    }
    $n=-1 unless $n;
    $avg = $sum/$n;
    $std = sqrt($sum2/$n-$avg*$avg);
    print O "# $avg ± $std\n";
    ($max_y,$max_x)=(-1,0);
    for my $k ($avg-$std .. $avg+$std) {	# 68.27%
	    next unless $v=$InsD{$k};
	    if ($max_y < $v) {
		    $max_x = $k;
		    $max_y = $v;
	    }
    }
    my $cutoff = $max_y / 1000;
    $cutoff = 3 if $cutoff<3;
    my ($diff, $Lc, $Rc);
    for my $k (keys %InsD) {
	    $v=$InsD{$k};
	    next if $v < $cutoff;
	    $diff = $k - $max_x;
	    if ($diff < 0) {
		    $Lsd += $v * $diff * $diff;
		    $Lc += $v;
	    }
	    elsif ($diff > 0) {
		    $Rsd += $v * $diff * $diff;
		    $Rc += $v;
	    }
    }
    $Lc=-1 unless $Lc;
    $Rc=-1 unless $Rc;
    $Lsd = sqrt($Lsd/$Lc);
    $Rsd = sqrt($Rsd/$Rc);
    print O "# +$Lsd -$Rsd\n#InsertSize\tCount\n";
    for my $k (sort {$a <=> $b} keys %InsD) {
	    $v=$InsD{$k};
	    print O "$k\t$v\n";
    }
    close O;
}

open NFO,'>',"$opt_o.info" or die "[x]Error opening $opt_o.info: $!\n";
if ($withPE) {
	my $p=sprintf "%.2f",100*$max_y/$n;
	$avg=int($avg*10+.5)/10;
	$Lsd=int($Lsd*100+.5)/100;
	$Rsd=int($Rsd*100+.5)/100;
	print NFO "#fmtS\tTotalReads\ttoPaired\ttoSingled\tModeIns(p%),Lsd,Rsd,InsAvg,STD\n";
	print NFO "Summary\t",join("\t",$InReads,$mapPair,$mapSingle,"$max_x($p %),$Lsd,$Rsd,$avg,$std"),"\n";
} else {
	print NFO "#fmtS\tTotalReads\tAlignment\n";
	print NFO "Summary\t",join("\t",$InReads,$mapSingle),"\n";
}
print NFO "\n#fmtC\tReadsOut\tBPOut\tMisSum\tTrimedReads\tTrimedBP\tmisMatchReads\tReads\@Hit\tBP\@Hit\tIndelReads\tBadLines\n";
print NFO join("\t",'ALL',$DatSum{'ReadsOut'},$DatSum{'BPOut'},$DatSum{'MisSum'},$DatSum{'TrimedReads'},$DatSum{'TrimedBP'},
	${&combineC($DatSum{'misMatch'})},${&combineJ($DatSum{'Hit9r'})},${&combineJ($DatSum{'Hit9bp'})},${&combineJ($DatSum{'Indel'})},$DatSum{'BadLines'}),"\n\n";# if exists $DatSum{'Hit9r'};
close NFO;

__END__
