#!/usr/bin/perl -w
use strict;
my ($file1, $file2) = @ARGV;
open IN,'<',"$file1" or die "$!";
my %hash;
while(<IN>){
    chomp;
    my($chr, $pos, $vaf) = (split "\t", $_)[0, 1, 2];
    $chr =~ s/chr//;
    $hash{$chr}{$pos} = $vaf;
}
close IN;

open IN,'<',"$file2" or die "$!";
while(<IN>){
    chomp;
    if($_ =~ m/^PeptidePostion/){
        print "$_\tDNAVAF\n";
        next;
    }
    my($chr, $pos) = (split "\t", $_)[9, 10];
    my @p = split /;/,$pos;
    my $vaf = 'na';
    foreach my $i (@p){
        $vaf = $hash{$chr}{$i} if(exists $hash{$chr}{$i});
    }
    print "$_\t$vaf\n";
}
