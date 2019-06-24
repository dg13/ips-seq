#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use Common qw(open_for_read2 open_for_write2 verbose checkFileExists checkDirExists);

my $usage = "$0 - Code to filter call set by FDR
    usage:
       $0 [OPTIONS] <inFile|stdin> <fdrThreshFile>
    options:
       -verbose=XXX - set verbosity level
       -XXX\n";

## Read options from CL
my $fdrThresh = 5;
GetOptions("verbose=i" => \$Common::verbosity, 'fdr=f' => \$fdrThresh);
chdir("/Users/dg13/Projects/ipsExome");

## Check CL args
&errAbort() unless @ARGV > 0;
my $inFile = $ARGV[0];
my $threshFile = $ARGV[1];
&checkFileExists($inFile) unless $inFile eq 'stdin';

## Read input files
my %fdr = &readFdrThresholds();
my $pval = $fdr{$fdrThresh}{'pval'};
my $prevLine;   
my $prevPval;
my $flag = 0;
my $pvalColNumber = 4;

*IN = &open_for_read2($inFile);
while(<IN>) {
    my @l = split;
    my $lastCol = scalar(@l) - 1;
    if($l[$lastCol]>2) { ## Check last column - format changes from older to newer CSQ version
	print STDERR "$0: Found > 2 nucleotide run on line $.\n";
	exit(1);
    }
    if($l[$lastCol]>1) {
	if($flag>0) {
	    if($l[$pvalColNumber]<$pval && $prevPval<$pval) { ## Pvalue for both bases of a dinuc must be below the threshold
		print $prevLine;
		print $_;
	    }
	    $flag = 0;
	} else {
	    $prevLine = $_;
	    $prevPval = $l[$pvalColNumber];
	    $flag++;
	}
    } else {
	print if $l[$pvalColNumber]<$pval;
    }
    
}
close(IN) unless $inFile eq 'stdin'; ## Avoids perl warning if we try and close STDIN;

## Functions

sub readFdrThresholds {
    my %dat;
    *IN = &open_for_read2($threshFile);
    while(<IN>) {
	next if /^fdr/;
	my @l = split;
	$dat{$l[0]}{'pval'} = $l[1];
    }
    close(IN);
    unless(exists $dat{$fdrThresh}) {
	print STDERR "$0: Specified FDR threshold ($fdrThresh) not found in $threshFile";
    }
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
