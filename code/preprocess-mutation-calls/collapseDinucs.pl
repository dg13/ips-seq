#!/usr/bin/perl -w

use lib "$ENV{'HOME'}/Projects/ipsExome/code";
use strict;
use Getopt::Long;
use Common qw(open_for_read2 open_for_write2 verbose checkFileExists checkDirExists);
use newCalls::newCalls;

my $usage = "$0 - Code to collapse dinucleotides to single bases, based on selection of most deleterious consequence
    usage:
       $0 [OPTIONS] <inFile|stdin>
    options:
       -verbose=XXX - set verbosity level
       -XXX\n";

## Read options from CL
GetOptions("verbose=i" => \$Common::verbosity);
chdir("/Users/dg13/Projects/ipsExome");

## Check CL args
&errAbort() unless @ARGV > 0;
my $inFile = $ARGV[0];
&checkFileExists($inFile) unless $inFile eq 'stdin';
my $csqFieldNo = 18;

## Read input files
my $n = 0;
my @l1;
*IN = &open_for_read2($inFile);
while(<IN>) {
    my @l = split;
    my $lastCol = scalar(@l) - 1;
    if($l[1] == 1110931) {
	print "";
    }
    if($.==44) {
	print "";
    }
    if($l[$lastCol]==2) {
	if($l[4]==179728586) {
	    print "";
	}
	$n++;
	if($n==2) {
	    my (@lens,@f1,@f2,@c1,@c2);
	    newCalls::getFields($l1[$csqFieldNo],\@f1,\@lens);
	    newCalls::getFields($l[$csqFieldNo],\@f2,\@lens);
	    newCalls::assignCsqCodes(\@f1,\@c1);
	    newCalls::assignCsqCodes(\@f2,\@c2);
	    if($c1[0][1] >= $c2[0][1]) { ## First base of dinuc has more deleterious consequences
##		printf "%s.1\n",join("\t",@l1);
	    } else {
##		printf "%s.2\n",join("\t",@l);
	    }
	    $n = 0;
	    @l1 = ();
	} else {
	    push(@l1,@l);
	}
    } else {
	## print;
    }
}
close(IN) unless $inFile eq 'stdin'; ## Avoids perl warning if we try and close STDIN;

## Functions 

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
