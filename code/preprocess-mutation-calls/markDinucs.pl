#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use Common qw(open_for_read2 open_for_write2 verbose checkFileExists checkDirExists);

my $usage = "$0 - Code to add column marking dinucleotides
    usage:
       $0 [OPTIONS] <inFile|stdin> <rlFile>
    options:
       -verbose=XXX - set verbosity level
       -XXX\n";

## Read options from CL
GetOptions("verbose=i" => \$Common::verbosity);
chdir("/Users/dg13/Projects/ipsExome");

## Check CL args
&errAbort() unless @ARGV > 1;
my $inFile = $ARGV[0];
my $rlFile = $ARGV[1];
&checkFileExists($inFile) unless $inFile eq 'stdin';

## Read input files
my $flag = 1;
my $timer = 0;
my %rl = &readRlFile($rlFile);
*IN = &open_for_read2($inFile);
while(<IN>) {
    chomp;
    my @l = split;
    if(exists $rl{$.}) {
	$flag = $rl{$.};
	$timer = $rl{$.};
	# if($flag==4) {
	#     print STDERR "$0: Error - quad nucleotide run found - wat do?\n";
	#     exit 1;
	# }
    }
    ## print "$_\t$flag\n" if $flag ==4;
    print "$_\t$flag\n";
    $timer--;
    $flag = 1 if $timer < 1;
}
close(IN) unless $inFile eq 'stdin'; ## Avoids perl warning if we try and close STDIN;

## Functions

sub readRlFile {
    my %dat;
    my $maxRun = `gunzip -c $rlFile | awk '\$3==1' | sort -k 2,2n | tail -n 1 | awk '{ print \$2 }'`;
    if($maxRun eq '') { ## No k-nuc runs found in file
	$dat{-1}++; ## Set to impossible value so not empty hash
	return %dat;
    }
    chomp($maxRun);
    if($maxRun > 4) {
	print STDERR "$0: Error - max number of sequential mutations is 2";
	exit(1);
    }
    *IN = &open_for_read2($rlFile);
    while(<IN>) {
	my @l = split;
	if($l[2] == 1) { ## There is a run of mutations separated by a single base
	    my $pos = $l[0] - ($l[1] - 1);
	    my $len = $l[1] + 1;
	    $dat{$pos} = $len;
	}
    }
    close(IN);
    return %dat;
}

sub usage {
    print STDERR $usage;
}

sub errAbort {
    ## Optional args:
    ## 1: Error message
    ## 2: 'q' suppresses usage message
    my $arg = shift;
    if(defined $arg) {
	print STDERR "Error: $arg\n";
    }
    undef $arg;
    $arg = shift(@_);
    if(defined $arg) {
	if($arg ne 'q') {
	    &usage();
	}
    } else {
	&usage();	
    }
    exit(1);
}
