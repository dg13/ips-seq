#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use Common qw(open_for_read2 open_for_write2 verbose checkFileExists checkDirExists);

my $usage = "$0 - Code to identify shared mutations between pairs of lines, where the site meets some p-value threshold in at least one of the sites
    Warning - to work, input file must contain just mutated sites from a pair of lines from the same donor, and must be sorted by chromosome then by position
    usage:
       $0 [OPTIONS] <inFile|stdin>
    options:
       -pval=XXX - pvalue threshold to define significance
       -verbose=XXX - set verbosity level
       -XXX\n";

## Read options from CL
chdir("/Users/dg13/Projects/ipsExome/");
my $fdrThresh = 0.05;
GetOptions("verbose=i" => \$Common::verbosity, 'pval=f' => \$fdrThresh);

## Check CL args
&errAbort() unless @ARGV > 0;
my $inFile = $ARGV[0];
&checkFileExists($inFile) unless $inFile eq 'stdin';

## Read input files
my @prevLine;
my %lineNames;
my $headerPrinted=0;
*IN = &open_for_read2($inFile);
while(<IN>) {
    my @l = split;
    $lineNames{$l[3]}++;
    if(!@prevLine) {
	push(@prevLine,@l);
	next;
    }
    if(!&linesIncompatible(\@l,\@prevLine)) { ## Tests whether the same position is tested in the correct pair of donors - some sites are filtered from one donor but not the other
	unless($headerPrinted) {
	    my @k = keys %lineNames;
	    print "#$k[0]\t$k[1]\n";
	    print "#chr\tpos\tref\talt\tpval1\tpval2\tfref\tfalt\tiref1\tialt1\tiref2\tialt2\n";
	    $headerPrinted++;
	}
	if(&eitherLineSignificant(\@l,\@prevLine)) {
	    &printLines(\@prevLine,\@l);
	} else {
	    print "LINE SPECIFIC $prevLine[9] $prevLine[10] $l[9] $l[10] $prevLine[0] $prevLine[1]\n";
	}
	@prevLine = ();
    } else {
	@prevLine = ();
	push(@prevLine,@l);
    }
}
close(IN) unless $inFile eq 'stdin'; ## Avoids perl warning if we try and close STDIN;

## Functions

sub printLines {
    my ($l1,$l2) = @_;
    print "$$l1[0]\t$$l1[1]\t$$l1[5]\t$$l1[6]\t$$l1[4]\t$$l2[4]\t$$l1[7]\t$$l1[8]\t$$l1[9]\t$$l1[10]\t$$l2[9]\t$$l2[10]\n";
}

sub eitherLineSignificant {
    my ($l1,$l2) = @_;
    if(($$l1[4] ne '.' && $$l1[4] < $fdrThresh) || ($$l2[4] ne '.' && $$l2[4] < $fdrThresh)) {
	return 1;
    }
    return 0;
}

sub linesIncompatible {
    my ($l1,$l2) = @_;
    my $ret = 0;
    $ret++ if $$l1[0] ne $$l2[0]; ## Diff chrs
    $ret++ if $$l1[1] != $$l2[1]; ## Diff positions
    $ret++ if $$l1[3] eq $$l2[3]; ## Same lines
    return $ret;
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
