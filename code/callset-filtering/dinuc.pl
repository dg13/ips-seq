#!/usr/bin/env perl
#
# Author: petr.danecek@sanger
#
# Filter dinucleotides. The input file format is described below.
#



use strict;
use warnings;
use Carp;

my $opts = parse_params();
if ( $$opts{allele_freqs} )
{
    correlate_dinuc_afs($opts);
}
elsif ( $$opts{positions} )
{
    print_dinuc_positions($opts);
}
else
{
    filter_dinuc($opts);
}

exit;

#--------------------------------

sub error
{
    my (@msg) = @_;
    if ( scalar @msg ) { confess @msg; }
    print 
        "Usage: dinuc [OPTIONS]\n",
        "Options:\n",
        "   -a, --allele-freqs            Correlation of allele frequencies\n",
        "   -p, --positions               Create a list of positions, suitable for check-dinuc-phase\n",
        "   -t, --threshold <float>       Fisher-test threshold\n",
        "   -h, -?, --help                This help message.\n",
        "\n";
    exit -1;
}
sub parse_params
{
    my $opts = {};
    while (defined(my $arg=shift(@ARGV)))
    {
        if ( $arg eq '-p' || $arg eq '--positions' ) { $$opts{positions} = 1; next; }
        if ( $arg eq '-a' || $arg eq '--allele-freqs' ) { $$opts{allele_freqs} = 1; next; }
        if ( $arg eq '-t' || $arg eq '--threshold' ) { $$opts{ft_th} = shift(@ARGV); next; }
        if ( $arg eq '-?' || $arg eq '-h' || $arg eq '--help' ) { error(); }
        error("Unknown parameter \"$arg\". Run -h for help.\n");
    }
    return $opts;
}


#  [1] chrom
#  [2] pos
#  [3] fibro
#  [4] ips
#  [5] Fisher test p-value
#  [6] REF
#  [7] ALT
#  [8] fibro: nREF
#  [9] fibro: nALT
#  [10] ips:   nREF
#  [11] ips:   nALT
sub correlate_dinuc_afs
{
    my ($opts) = @_;
    my %prev = ();
    while (my $line=<STDIN>)
    {
        if ( $line=~/^#/ ) { next; }
        my @vals = split(/\t/,$line);
        my $mut = "$vals[5]$vals[6]";

        # Test if things change much if only UV-consistent dinucs are included
        #   if ( $mut ne 'GA' && $mut ne "CT" ) { next; }

        my $ft  = $vals[4];
        if ( exists($$opts{ft_th}) && $$opts{ft_th} < $ft ) { next; }
        my $ips = $vals[3];
        my $fib = $vals[2];
        my $pos = $vals[1];
        my $chr = $vals[0];
        my $af_fib = $vals[8]/($vals[7]+$vals[8]);
        my $af_ips = $vals[10]/($vals[9]+$vals[10]);
        my $key = "$fib.$ips";
        if ( exists($prev{$key}) )
        {
            if ( $prev{$key}{chr} ne $chr ) { goto done; }
            if ( $prev{$key}{pos}+1 != $pos ) { goto done; }
            my $max_fib = $af_fib > $prev{$key}{af_fib} ? $af_fib : $prev{$key}{af_fib};
            my $max_ips = $af_ips > $prev{$key}{af_ips} ? $af_ips : $prev{$key}{af_ips};
            my ($af1,$af2);
            if ( $max_fib > $max_ips ) { $af1 = $af_ips; $af2 = $prev{$key}{af_ips}; }
            else { $af1 = $af_fib; $af2 = $prev{$key}{af_fib}; }
            if ( $af1 or $af2 ) { print "$af1\t$af2\n"; }
        }
done:
        $prev{$key}{mut} = $mut;
        $prev{$key}{chr} = $chr;
        $prev{$key}{pos} = $pos;
        $prev{$key}{ft}  = $ft;
        $prev{$key}{rec} = $line;
        $prev{$key}{af_fib} = $af_fib;
        $prev{$key}{af_ips} = $af_ips;
    }
}

#  [1] chrom
#  [2] pos
#  [3] fibro
#  [4] ips
#  [5] Fisher test p-value
#  [6] REF
#  [7] ALT
#  [8] fibro: nREF
#  [9] fibro: nALT
#  [10] ips:   nREF
#  [11] ips:   nALT
sub print_dinuc_positions
{
    my ($opts) = @_;
    my %prev = ();
    while (my $line=<STDIN>)
    {
        if ( $line=~/^#/ ) { next; }
        my @vals = split(/\t/,$line);
        my $alt = $vals[6];
        my $mut = "$vals[5]$vals[6]";

        my $ft  = $vals[4];
        if ( exists($$opts{ft_th}) && $$opts{ft_th} < $ft ) { next; }
        my $ips = $vals[3];
        my $fib = $vals[2];
        my $pos = $vals[1];
        my $chr = $vals[0];
        my $af_fib = $vals[8]/($vals[7]+$vals[8]);
        my $af_ips = $vals[10]/($vals[9]+$vals[10]);
        my $key = "$fib.$ips";
        if ( exists($prev{$key}) )
        {
            if ( $prev{$key}{chr} ne $chr ) { goto done; }
            if ( $prev{$key}{pos}+1 != $pos ) { goto done; }
            my $max_fib = $af_fib > $prev{$key}{af_fib} ? $af_fib : $prev{$key}{af_fib};
            my $max_ips = $af_ips > $prev{$key}{af_ips} ? $af_ips : $prev{$key}{af_ips};
            my ($af1,$af2);
            if ( $max_fib > $max_ips ) { print "$fib\t$chr\t$prev{$key}{pos}\t$alt\n"; }
            else { print "$ips\t$chr\t$prev{$key}{pos}\t$alt\n"; }
        }
done:
        $prev{$key}{alt} = $alt;
        $prev{$key}{mut} = $mut;
        $prev{$key}{chr} = $chr;
        $prev{$key}{pos} = $pos;
        $prev{$key}{ft}  = $ft;
        $prev{$key}{rec} = $line;
        $prev{$key}{af_fib} = $af_fib;
        $prev{$key}{af_ips} = $af_ips;
    }
}

#  [1] chrom
#  [2] pos
#  [3] fibro
#  [4] ips
#  [5] Fisher test p-value
#  [6] REF
#  [7] ALT
sub filter_dinuc
{
    my ($opts) = @_;
    my %prev = ();
    while (my $line=<STDIN>)
    {
        if ( $line=~/^#/ ) { print $line; next; }
        my @vals = split(/\t/,$line);
        my $mut = "$vals[5]$vals[6]";
        if ( $mut ne 'GA' && $mut ne "CT" ) { next; }
        my $ft  = $vals[4];
        my $ips = $vals[3];
        my $fib = $vals[2];
        my $pos = $vals[1];
        my $chr = $vals[0];
        my $key = "$fib.$ips";
        if ( exists($prev{$key}) )
        {
            if ( $prev{$key}{mut} ne $mut ) { goto done; }
            if ( $prev{$key}{chr} ne $chr ) { goto done; }
            if ( $prev{$key}{pos}+1 != $pos ) { goto done; }
            if ( $prev{$key}{ft} > $ft ) { $prev{$key}{ft} = $ft; }
            if ( exists($$opts{ft_th}) && $$opts{ft_th} < $prev{$key}{ft} ) { next; }
            print $prev{$key}{rec};
            print $line;
            delete($prev{$key});
            next;
        }
done:
        $prev{$key}{mut} = $mut;
        $prev{$key}{chr} = $chr;
        $prev{$key}{pos} = $pos;
        $prev{$key}{ft}  = $ft;
        $prev{$key}{rec} = $line;
    }
}


