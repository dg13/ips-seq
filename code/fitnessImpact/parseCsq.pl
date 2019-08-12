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

chdir("/Users/dg13/Projects/ipsExome");

## Read options from CL
my $print_alleles = 0;
GetOptions("verbose=i" => \$Common::verbosity, "print-alleles" => \$print_alleles);

## Check CL args
&errAbort() unless @ARGV > 0;
my $inFile = $ARGV[0];
&checkFileExists($inFile) unless $inFile eq 'stdin';
my $sep = "\t";

## Read input files
*IN = &open_for_read2($inFile);
while(<IN>) {
    my @l = split;
    if(scalar(@l)!=6 && $print_alleles==0) {
	my $n = scalar(@l);
	printf STDERR "Error unexpected line length ($n fields)\n";
	exit 1;
    }
    if($l[1]==13418) {
	print "";
    }
    my $ref = my $alt = "";
    if($print_alleles>0) {
	if(scalar(@l)!=8) {
	    my $n = scalar(@l);
	    printf STDERR "Error unexpected line length ($n fields)\n";
	    exit 1;
	}
	$alt = pop(@l);
	$ref = pop(@l);
    }
    my $bcsq = pop(@l);
    my $csq = pop(@l);
    my $siftCount = my $pphenCount = 0;
    my $sift = my $pphen = 0;
    my (%bcsqGene,%bcsqImpact,%csqGene,%csqImpact,%bcsq);
    &processCsq($csq,\%csqGene,\%csqImpact,\$siftCount,\$pphenCount,\$sift,\$pphen);
    &processBcsq($bcsq,\%bcsqGene,\%bcsqImpact);
    ## next unless $pphenCount > 0 || $siftCount > 0; ## Exclude lines without any impact annotation
    my $out = join($sep,@l);
    my $csqGene = join(",",keys %csqGene);
    my $csqImpact = join(",",keys %csqImpact);
    my $bcsqGene = join(",",keys %bcsqGene);
    my $bcsqImpact = join(",",keys %bcsqImpact);
    print "$out$sep$csqGene$sep$csqImpact$sep$bcsqGene$sep$bcsqImpact$sep$sift$sep$pphen\t$ref\t$alt\n";
}
close(IN) unless $inFile eq 'stdin'; ## Avoids perl warning if we try and close STDIN;

## Functions

sub processBcsq {
    my ($csq,$gene,$impact) = @_;
    if($csq eq "." || $csq =~ /^@[0-9]*$/) {
	$$impact{'NA'}++;
	$$gene{'NA'}++;
    } else {
	my @csq = split(",",$csq);
	foreach my $f (@csq) {
	    next if $f =~ /^@/;
	    my @f = split("\\|",$f,-1);
	    $$impact{$f[0]}++;
	    $$gene{$f[1]}++;
	}
    }
    unless(%{ $gene }) {
	$$gene{'NA'}++;
    }
    unless(%{ $impact }) {
	$$impact{'NA'}++;
    }
    return 0;
}

sub processCsq {
    my ($csq,$gene,$impact,$siftCount,$pphenCount,$sift,$pphen) = @_;
    my @csq = split(",",$csq);
    my $flag=0;
    foreach my $f (@csq) {
	next if $f eq '.';
	my @f = split("\\|",$f,-1);
	if(scalar(@f)!=24) {
	    print STDERR $_;
	    print STDERR "Syntax error\n";
	    exit;
	}
	if($f[22] ne '') {
	    $flag++;
	    my @n = split("[()]",$f[22]);
	    $$sift += $n[1];
	    $$siftCount++;
	}
	if($f[23] ne '') {
	    $flag++;
	    my @n = split("[()]",$f[23]);
	    $$pphen += $n[1];
	    $$pphenCount++;
	}
	my @impact = split("&",$f[1]);
	foreach my $imp (@impact) {
	    $$impact{$imp}++;
	}
	$f[3] = 'NA' if $f[3] eq '';
	$$gene{$f[3]}++;
    }
    ## Get average across all annotated impacts
    if($$siftCount > 0) {
	$$sift /= $$siftCount;
    } else {
	$$sift = 'NA';
    }
    if($$pphenCount > 0) {
	$$pphen /= $$pphenCount;
    } else {
	$$pphen = 'NA';
    }
    unless(%{ $gene }) {
	$$gene{'NA'}++;
    }
    unless(%{ $impact }) {
	$$impact{'NA'}++;
    }
    return 0;
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
