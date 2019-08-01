#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use List::MoreUtils qw{uniq};

my ($in, $out, $popmap, $cov, $help, $pops, $cnt, $error, $type, %map_ind);
my ($in_fh);
GetOptions (
    "in=s"      =>\$in,
    "out=s"     =>\$out,
    "popmap=s"  =>\$popmap,
    "cov=f"     =>\$cov,
    "error=f"   =>\$error,
    "type=s"    =>\$type,
    "help:1"    =>\$help
) or die "$!";

my $usage = '
    -i: input GT file.
    -o: output file.
    -p: popmap file.
    -c: coverage.
    -e: error rate. 
    -h: help message.
';
$cov   = defined $cov ? $cov : 0.6;
$error = defined $error ? $error : 0.05;
die $usage if ($help || !($in && $out && $popmap));
if (substr($in, -2) eq 'gz') {
    open( $in_fh, "gzip -dc $in|") or die "$!";
} else {
    open($in_fh, "$in") or die "$!";
}

($pops, $cnt) = parse_popmap($popmap);

print STDERR $cnt->{'male'}+$cnt->{'female'}, " individuals, $cnt->{'male'} males, $cnt->{'female'} females\n";
open(my $out_fh, ">$out") or die "$!";
print $out_fh join("\t", 'chrom', 'pos', 'female_'.$cnt->{'female'}, 'male_'.$cnt->{'male'}, 'het_f', 'het_m');

my $th_male   = $cov   >= 1 ? $cov   : $cov*$cnt->{'male'};
my $e_male    = $error >= 1 ? $error : $error*$cnt->{'female'};
my $th_female = $cov   >= 1 ? $cov   : $cov*$cnt->{'female'};
my $e_female  = $error >= 1 ? $error : $error*$cnt->{'male'};

my $vcf_header = 0;
my $sex_loci_cnt = 0;
while (<$in_fh>) {
    
    # GT from vcftools
    # CHROM POS IND...
    my ($chrom, $pos, $id, $i, @females, @males, @parts, @GT, %het, %hom, %all, $hom_gt);
    next if /^#|^$/;
    chomp;
    
    if ($vcf_header == 0) {
        $vcf_header = 1;
        $i = 0;
        @parts  = split;
        map {$map_ind{++$i} = $_;print $out_fh "\t$_";} @parts[2..$#parts];
        print $out_fh "\n";
        next;
    }
    
    @parts = split;
    $chrom = $parts[0];
    $pos   = $parts[1];
    $i     = 0;
    
    foreach(@parts[2..$#parts]) {
        push @GT, $_;
        ++$i;
        my $gt  = $_;
        next if $gt eq './.';
        my $ind = $map_ind{$i};
        my $sex = $pops->{$ind};
        if ($gt eq '0/1') {
            $het{$sex}++;
        } else {
            $hom_gt->{$sex}->{$gt}++; # gt type for hom.
            $hom{$sex}++;
        }
        $all{$sex}++;
    }
    
    my $f     = $het{'female'} ? $het{'female'} : 0;
    my $m     = $het{'male'} ? $het{'male'} : 0;
    
    if (($f >= $th_female && $m <= $e_female && scalar(keys %{$hom_gt->{'male'}}) == 1) || 
    ($m >= $th_male && $f <= $e_male && scalar(keys %{$hom_gt->{'female'}}) == 1)) {
        # putative sex loci.
        my $f_het = sprintf "%.5f", $f/$all{'female'};
        my $m_het = sprintf "%.5f", $m/$all{'male'};
        print $out_fh join("\t", $chrom, $pos, $f, $m, $f_het, $m_het, @GT), "\n";
        $sex_loci_cnt++;
    }
    
}
close $in_fh;
close $out_fh;
print STDERR "Total $sex_loci_cnt putative sex specific loci.\n";

sub parse_popmap {
    
    my $pop = shift;
    my ($pops, $cnt);
    open(my $in_fh, $pop) or die "$!";
    while(<$in_fh>) {
        next if /^#|^$/;
        $_ =~ s/\r//g;
        chomp;
        # indv sex
        my @parts = split;
        die if scalar(@parts) ne 2;
        my $ind   = $parts[0];
        my $sex   = lc($parts[1]);
        $pops->{$ind} = $sex;
        $cnt->{$sex}++;
    }
    close $in_fh;
    return $pops, $cnt;
    
}

sub max_cnt {

    my $hash = shift;
    return 0 if not %$hash;
    my @cnt = sort {$b<=>$a} values %$hash;
    return $cnt[0];
}
