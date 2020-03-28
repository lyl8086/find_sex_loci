#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use List::MoreUtils qw{uniq};

my ($in, $out, $popmap, $cov, $help, $pops, $cnt, $error, $type);
my ($in_fh, $loci, $v);
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
    -i: input catalog file.
    -o: output file.
    -p: popmap file.
    -c: coverage.
    -e: error rate.
    -t: haplotype file. 
    -h: help message.
';
$cov   = defined $cov ? $cov : 0.8;
$error = defined $error ? $error : 0.05;
die $usage if ($help || !($in && $out && $popmap));
if (substr($in, -2) eq 'gz') {
    open( $in_fh, "gzip -dc $in|") or die "$!";
} else {
    open($in_fh, "$in") or die "$!";
}

($pops, $cnt) = parse_popmap($popmap);

print STDERR $cnt->{'male'}+$cnt->{'female'}," individuals, ", "$cnt->{'male'} males, $cnt->{'female'} females\n";

open(my $out_fh, ">$out.catalog") or die "$!";
print $out_fh join("\t", 'id', 'female_'.$cnt->{'female'}, 'male_'.$cnt->{'male'}, 'seq'), "\n";

my $th_male   = $cov   >= 1 ? $cov   : $cov*$cnt->{'male'};
my $e_male    = $error >= 1 ? $error : $error*$cnt->{'female'};
my $th_female = $cov   >= 1 ? $cov   : $cov*$cnt->{'female'};
my $e_female  = $error >= 1 ? $error : $error*$cnt->{'male'};
while (<$in_fh>) {
    
    # Cstacks
    # 0       12      1               0       +       consensus       0       23_797,37_75507,38_127628,41_42213,46_133241,54_53470,55_75236,56_147331,68_67432,73_39739,70_137313,71_75264,78_120962,81_27836,82_24828,135_138769,137_31077,139_90684,50_18824,49_15174,51_26070,53_26528,57_105955,61_61442,58_49263,59_85479,60_36774,62_105377,63_131753,64_65132,65_60678,67_69778,69_28639,72_65659,74_47377 AATTCAGCTTCTGTTTTTTTTTTTATTATTATAAGAAATATGACTCAGCCGAAAGCTTTAAACTGTATGAAATGTTAATAACTAATCCCAAAATCCAAACTTCTGTTCCACCCCACAGAGGTGTC        0       0       0       0
    # sql_id sample_id  locus_id    chromesome  base    strand  model_type    stacks_components   seq_id  seq   Deleveraged_Flag    Blacklisted_Flag    Lumberjackstack_Flag    likelihood
    #
    my (@males, @females, $id, $seq, $inds);
    next if /^#|^$/;
    chomp;
    my @parts  = split(/\t/);
    if (!$v) {
        # only the use the first line.
        my $npart = @parts;
        if ($npart == 14) {$v = 1;}
        elsif ($npart == 9) {$v = 2;}
        else {die "Stacks tags files error!";}
    }
   
    if ($v == 1) {
        $id   = $parts[2];
        $seq  = $parts[9];
        $inds = $parts[8];
    } elsif ($v == 2) {
        $id   = $parts[1];
        $seq  = $parts[5];
        $inds = $parts[4];
    } else {
        die "Stacks tags files error!";
    } 
    $loci->{$id} = { 'id'=>$id, 'seq'=>$seq};
    my @indivs = split(/,/, $inds);
    foreach(@indivs) {
        my $ind = (split(/_/, $_))[0];
        my $sex = $pops->{$ind}->{'sex'};
        push @females, $ind if $sex eq 'female';
        push @males, $ind if $sex eq 'male';
    }
    @females = uniq(@females);
    @males   = uniq(@males);
    my $f = @females ? @females : 0;
    my $m = @males ? @males : 0;
    
    if (($f >= $th_female && $m <= $e_female) || ($m >= $th_male && $f <= $e_male)) {
        # putative sex loci.
        my @fake = $f > $m ? @males : @females;
        print $out_fh join("\t", $id, $f, $m, $seq, map{$pops->{$_}->{'ind'};} @fake), "\n";
    }
    
}
close $in_fh;
close $out_fh;

# haplotype file.
if (defined $type) {
open($in_fh, "$type") or die "$!";
open($out_fh, ">$out.hap") or die "$!";
open(my $out_fh_2, ">$out.hap.het") or die "$!";
my $tmp = <$in_fh>;
my $grp;
$tmp =~ s/[\r\n]//g;
my @head = split(/\t/, $tmp);

foreach (@head[2..$#head]) {
    my $sex = $pops->{'sex'}->{$_};
    push @{$grp->{$sex}}, $_;
    
}

print join("\t", 'cat_id', 'female', 'male', @{$grp->{'female'}}, '|', @{$grp->{'male'}}), "\n";
while(<$in_fh>) {
    chomp;
    my @parts   = split(/\t/);
    my $cat_id  = $parts[0];
    my $hap_cnt = $parts[1]; next if $hap_cnt == 1; # 
    my $i       = 2;
    my (%het, %geno, %male_hap, %female_hap);
    my (@males, @females);
    foreach my $hap (@parts[2..$#parts]) {
        
        my $indiv = $head[$i];
        my $sex   = $pops->{'sex'}->{$indiv};
        $het{$sex}++ if $hap =~ /\//;
        $geno{$indiv} = $hap;
        $i++;
        next if $hap eq '-';
        if ($sex eq 'male') {
            push @males, $indiv;
            $male_hap{$hap}++;
        } elsif ($sex eq 'female') {
            push @females, $indiv;
            $female_hap{$hap}++;
        } else {
            die "error in header...\n";
        }
    }
    
    # f: num of females.
    # m: num of males.
    my $f = @females ? @females : 0;
    my $m = @males ? @males : 0;
    next if ($f + $m == 0);
    if (($f >= $th_female && $m <= $e_female) || ($m >= $th_male && $f <= $e_male)) {
        # exclusive.
        my @fake = $f > $m ? @males : @females;
        print $out_fh join("\t", $cat_id, $f, $m, $loci->{$cat_id}->{'seq'}, @fake), "\n";
    }
    
    $m = max_cnt(\%male_hap);
    $f = max_cnt(\%female_hap);
    
    next if ($f < $th_female || $m < $th_male); # not enough coverage.
    
    $f = $het{'female'} ? $het{'female'} : 0;
    $m = $het{'male'} ? $het{'male'} : 0;
    
    if (($f >= $th_female && $m <= $e_female) || ($m >= $th_male && $f <= $e_male)) {
        my @fake = $f > $m ? @males : @females;
        print $out_fh_2 join("\t", $cat_id, $f, $m, $loci->{$cat_id}->{'seq'}, map {$_.':'.$geno{$_}} @fake), "\n";
        print join("\t", $cat_id, $f, $m, map {$geno{$_}} @{$grp->{'female'}}), "\t"; 
        print join("\t", '|', map {$geno{$_}} @{$grp->{'male'}}), "\n";
    }
}
}

sub parse_popmap {
    
    my $pop = shift;
    my ($pops, $cnt);
    open(my $in_fh, $pop) or die "$!";
    while(<$in_fh>) {
        next if /^#|^$/;
        $_ =~ s/\r//g;
        chomp;
        # indv  id  sex
        # id: individual number in stacks.
        my @parts = split;
        die if scalar(@parts) ne 3;
        my $ind   = $parts[0];
        my $id    = $parts[1];
        my $sex   = lc($parts[2]);
        $pops->{'sex'}->{$ind} = $sex; # for type II.
        if ($sex eq 'female') {
            $pops->{$id}->{'sex'} = 'female';
            $pops->{$id}->{'ind'} = $ind;
            $cnt->{'female'}++;
        } elsif ($sex eq 'male') {
            $pops->{$id}->{'sex'} = 'male';
            $pops->{$id}->{'ind'} = $ind;
            $cnt->{'male'}++;
        }
        
        
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
