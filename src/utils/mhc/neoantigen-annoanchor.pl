#!/usr/bin/env perl 
use strict;
use warnings;
use utf8;

my $file = shift;
open IN, "<$file" or die "$!";
my %hash;
while(<IN>){
    chomp;
    next if($_ =~ m/Genotype/);
    my($fac, $length, $sites) = split "\t", $_;
    $hash{$fac}{$length} = $sites if($sites && $sites ne "");
}

my $annofile = shift;
open INA, "<$annofile" or die "$!";
while(<INA>){
    chomp;
    if($_ =~ m/^PeptidePostion/){
        print "$_\tAnchorFlag\tPsite\n";
        next;
    }
    my($pos, $HLA, $Peptide, $original) = (split "\t", $_)[0, 1, 2, 17];
    my $flag = 0;
    my $length = length($Peptide);
    my @orig = split /;/,$original;
    if(exists $hash{$HLA}{$length} && $hash{$HLA}{$length} ne ""){
        my @ps = split ",", $hash{$HLA}{$length};
        foreach my $pi(@ps){
            foreach my $pr (@orig){
                my $p = $pr - $pos + 1;
                $flag = $p if($p == $pi);
            }
        }
    }
    my $relpos = '';
    foreach my $i (@orig){
        my $p = $i - $pos + 1;
        if(length($relpos) == 0){
            $relpos = $p;
        } else{
            $relpos .= ";$p";
        }
    }
    print "$_\t$flag\t$relpos\n";
}


#===============================================================================
#
#         FILE: matrix2table.pl
#
#        USAGE: ./matrix2table.pl  
#
#  DESCRIPTION: 
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Jiaqian Wang, wangjq@yucebio.com
# ORGANIZATION: YuceBio
#      VERSION: 1.0
#      CREATED: 2015年11月11日 12时10分52秒
#     REVISION: ---
#===============================================================================


