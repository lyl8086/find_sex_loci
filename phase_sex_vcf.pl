#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use List::MoreUtils qw{uniq};

my ($in, $out, $popmap, $cov, $help, $pops, $cnt, $error, $mhet, %map_ind);
my ($in_fh);
GetOptions (
    "in=s"      =>\$in,
    "out=s"     =>\$out,
    "popmap=s"  =>\$popmap,
    "cov=f"     =>\$cov,
    "error=f"   =>\$error,
    "zhd=f"     =>\$mhet,
    "help:1"    =>\$help
) or die "$!";

my $usage = '
    -i: input vcf file.
    -o: output vcf file.
    -p: popmap file.
    -c: coverage.
    -e: error rate.
    -z: minimum heterozygosity.
    -h: help message.
';
$cov   = defined $cov ? $cov : 0.6;
$error = defined $error ? $error : 0.05;
$mhet  = defined $mhet ? $mhet : 0.8;
die $usage if ($help || !($in && $out && $popmap));
if (substr($in, -2) eq 'gz') {
    open( $in_fh, "gzip -dc $in|") or die "$!";
} else {
    open($in_fh, "$in") or die "$!";
}

($pops, $cnt) = parse_popmap($popmap);

print STDERR $cnt->{'male'}+$cnt->{'female'}, " individuals, $cnt->{'male'} males, $cnt->{'female'} females\n";
open(my $out_fh, ">$out") or die "$!";

my $th_male   = $cov   >= 1 ? $cov   : $cov*$cnt->{'male'};
my $e_male    = $error >= 1 ? $error : $error*$cnt->{'female'};
my $th_female = $cov   >= 1 ? $cov   : $cov*$cnt->{'female'};
my $e_female  = $error >= 1 ? $error : $error*$cnt->{'male'};

my $vcf_header = 0;
my $sex_loci_cnt = 0;
my @header;
while (<$in_fh>) {
    #####
    if (/^##/) {
        print $out_fh $_;
        next;
    }
    chomp;
    if(/^#/ && $vcf_header == 0) {
        @header = split(/\t/, $_);
        $vcf_header = 1;
        print $out_fh join("\t",@header[0..8], "het"),"\n";
    }
    my ($chrom, $pos, $id, $ref, $alt, $qual, $filter, $infoo, $format, @genos) = split(/\t/, $_);
    my ($i, $fmt, @females, @males, @parts, $GT, %het, %hom, %all, $hom_gt);
    next if /^#|^$/;
    my @formats = split(/:/, $format);              # Order of format:
    map { $fmt->{$formats[$_]} = $_;} 0..$#formats; # Geno => Order.
    foreach my $geno (@genos) {
        my @geno = split(/:/, $geno);
        ++$i;
        my $gt  = $geno[$fmt->{'GT'}];
        next if $gt eq './.';
        my $ind = $header[$i+8];
        my $sex = $pops->{$ind};
        if ($gt eq '0/1' || $gt eq '0|1' || $gt eq '1|0' || $gt eq '1/0') {
            $het{$sex}++;
        } else {
            $hom_gt->{$sex}->{$gt}++; # gt type for hom.
            $hom{$sex}++;
        }
        $all{$sex}++;
    }
    
    my $f     = $het{'female'} ? $het{'female'} : 0;
    my $m     = $het{'male'} ? $het{'male'} : 0;
    my @male_hom_gt    = (keys %{$hom_gt->{'male'}});
    my @female_hom_gt = (keys %{$hom_gt->{'female'}});
    my $f_het = sprintf "%.5f", $f/$all{'female'};
    my $m_het = sprintf "%.5f", $m/$all{'male'};
    my $info  = join(";","hetf:$f_het","hetm:$m_het","nf:$all{'female'}","nm:$all{'male'}");
    
    if ($f >= $th_female && $m <= $e_female && scalar(@male_hom_gt) == 1 && $f_het >= $mhet) {
        # ZW.
        my $ps = '01';
        if ($male_hom_gt[0] eq '0/0' || $male_hom_gt[0] eq '0|0') {
            $GT = '0|1';
        } else {
            $GT = '1|0';
        }
        print $out_fh join("\t", $chrom, $pos, $id, $ref, $alt, $qual, $filter, "ZW;$info", 'GT:PS', "$GT:$ps"), "\n";
        $sex_loci_cnt++;
    } 
    if ($m >= $th_male && $f <= $e_male && scalar(@female_hom_gt) == 1 && $m_het >= $mhet) {
        # XY.
        my $ps = '02';
        if ($female_hom_gt[0] eq '0/0' || $female_hom_gt[0] eq '0|0') {
            $GT = '0|1';
        } else {
            $GT = '1|0';
        }
        print $out_fh join("\t", $chrom, $pos, $id, $ref, $alt, $qual, $filter, "XY;$info", 'GT:PS', "$GT:$ps"), "\n";
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
