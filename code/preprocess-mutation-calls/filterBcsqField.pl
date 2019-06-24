#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use Common qw(open_for_read2 open_for_write2 verbose checkFileExists checkDirExists);

my $usage = "$0 - Code to pick the highest impact version of a mutation from multiple options in the BCSQ field
    usage:
       $0 [OPTIONS] <inFile|stdin>
    options:
       -verbose=XXX - set verbosity level
       -XXX\n";

## Read options from CL
GetOptions("verbose=i" => \$Common::verbosity);
## chdir("/Users/dg13/Projects/ipsExome");

## Check CL args
&errAbort() unless @ARGV > 0;
my $inFile = $ARGV[0];
&checkFileExists($inFile) unless $inFile eq 'stdin';
my $csqFieldNo = 18;

## Read input files
*IN = &open_for_read2($inFile);
while(<IN>) {
    my @l = split;
    $l[$csqFieldNo] = &processBcsqField($l[$csqFieldNo]);
    printf "%s\n",join("\t",@l);
}
close(IN) unless $inFile eq 'stdin'; ## Avoids perl warning if we try and close STDIN;

## Functions

sub processBcsqField {
    my $csq = shift;
    if($csq =~ /,/) {
	my (@fields,@lens,@codes);
	&getFields($csq,\@fields,\@lens);
	&assignCsqCodes(\@fields,\@codes);
	$csq = &getMostDeleteriousCsq(\@fields,\@codes);
    }
    return $csq;
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

sub getFields {
    my ($csq,$ss,$l) = @_;
    my @s = split(/,/,$csq);
    my $i = 0;
    foreach my $f (@s) {
	push(@{ $$ss[$i] },split(/\|/,$f));
	$$l[$i] = scalar(@{ $$ss[$i] });
	$i++;
    }
}

sub getCode {
    ## Finds a CSQ label, and assigns a numeric code based on order of deleteriousness
    my ($str) = @_;
    my $ret = 0;
    if($str =~ /&/) {
	my @s = split(/&/,$str);
	$str = $s[0];
    }
    $str =~ s/\*//;
    if($str eq '.') {
	$ret = 0;
    } elsif($str eq 'non_coding') {
	$ret = 0;
    } elsif($str =~ /^@[0-9]+/) {
	$ret = 0;
    } elsif($str eq 'intron') {
	$ret = 1;
    } elsif($str eq 'synonymous') {
	$ret = 3;
    } elsif($str eq 'stop_retained') {
	$ret = 2;
    } elsif($str eq '3_prime_utr') {
	$ret = 4;
    } elsif($str eq '5_prime_utr') {
	$ret = 5;
    } elsif($str eq 'splice_region') {
	$ret = 6;
    } elsif($str eq 'missense') {
	$ret = 7;
    } elsif($str eq 'splice_donor') {
	$ret = 8;
    } elsif($str eq 'splice_acceptor') {
	$ret = 9;
    } elsif($str eq 'start_lost') {
	$ret = 10;
    } elsif($str eq 'stop_lost') {
	$ret = 11;
    } elsif($str eq 'stop_gained') {
	$ret = 12;
    } else {
	print STDERR "$0: Unknown consequence $str";
	exit(1);
    }
}

sub assignCsqCodes {
    my ($fi,$co) = @_;
    my $i=0;
    foreach my $fi (@$fi) {
	my $code = &getCode($fi->[0]);
	$$co[$i][0] = $i;
	$$co[$i][1] = $code;
	$i++;
    }
}

sub getMostDeleteriousCsq {
    my ($f,$c) = @_;
    @$c = sort { $b->[1] <=> $a->[1] } @$c;
    return join("|",@{ $$f[$$c[0][0]] });
}
