#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use Common qw(open_for_read2 open_for_write2 verbose checkFileExists checkDirExists);

my $usage = "$0 - Code to count the number of donors a variable site is found in, for subsequent removal
    usage:
       $0 [OPTIONS] <inFile|stdin>
    options:
       -verbose=XXX - set verbosity level
       -XXX\n";

chdir("/Users/dg13/Projects/ipsExome");

## Read options from CL
GetOptions("verbose=i" => \$Common::verbosity);

## Check CL args
&errAbort() unless @ARGV > 0;
my $inFile = $ARGV[0];
&checkFileExists($inFile) unless $inFile eq 'stdin';
## Read input files
my %count;
*IN = &open_for_read2($inFile);
while(<IN>) {
    my @l = split;
    $count{$l[0]}{$l[1]}{$l[2]}++;
}
close(IN) unless $inFile eq 'stdin'; ## Avoids perl warning if we try and close STDIN;

foreach my $chr (sort keys %count) {
    foreach my $pos (sort { $a <=> $b } keys %{ $count{$chr} }) {
	my $n = scalar keys %{ $count{$chr}{$pos} };
	if($n>1) {
	    print "$chr\t$pos\t$n\n";
	}
    }
}

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
