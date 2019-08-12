#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use Common qw(open_for_read2 open_for_write2 verbose checkFileExists checkDirExists);

my $usage = "$0 - Code to parse and extract info from the CSQ field in petrs VCFs
    usage:
       $0 [OPTIONS] <inFile|stdin>
    options:
       -verbose=XXX - set verbosity level
       -XXX\n";

## Read options from CL
my $mafthresh = 0.01;
GetOptions('maf=f' => \$mafthresh);

## Check CL args
&errAbort() unless @ARGV > 0;
my $inFile = $ARGV[0];
&checkFileExists($inFile) unless $inFile eq 'stdin';
my $sep = "\t";
## Read input files
*IN = &open_for_read2($inFile);
while(<IN>) {
    my @l = split;
    my $af1kg = $l[2];
    my $afexac = $l[3];
    if($af1kg =~ /,/) {
	my @tmp = split(",",$af1kg);
	$af1kg = $tmp[0];
    }
    if($afexac =~ /,/) {
	my @tmp = split(",",$afexac);
	$afexac = $tmp[0];
    }
    next if $af1kg eq '.';
    next if $afexac eq '.';
    next if $af1kg <= $mafthresh || $afexac <= $mafthresh;
    print;
}
close(IN) unless $inFile eq 'stdin'; ## Avoids perl warning if we try and close STDIN;
