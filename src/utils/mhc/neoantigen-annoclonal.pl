#!/usr/bin/perl -w
use strict;
my ($file1, $file2, $i) = @ARGV;
my %hash;
if(-e $file1){
    open IN,'<',"$file1" or die "$!";
    while(<IN>){
        chomp;
        next if($_ =~ /^Chr/);
        my($chr, $pos, $vaf) = (split "\t", $_)[0, 1, 12];
        $chr =~ s/^chr//;
        $hash{$chr}{$pos} = $vaf;
    }
    close IN;
}

open IN,'<',"$file2" or die "$!";
while(<IN>){
    chomp;
    if($_ =~ m/^NeoRank/){
        print "$_\tCCF\tClonal\n";
        next;
    }
    my($chr, $poss) = (split "\t", $_)[13, 14];
    $chr =~ s/^chr//;
    my $pos = $poss;
    if($poss =~ /;/){
        $pos = (split /;/, $poss)[0];
    }
    my $clonal = "na";
    $clonal = "subclonal" if(exists $hash{$chr}{$pos} && $hash{$chr}{$pos} < 0.5);
    $clonal = "clonal" if(exists $hash{$chr}{$pos} && $hash{$chr}{$pos} >= 0.5);
    print "$_\t$hash{$chr}{$pos}\t$clonal\n" if(exists $hash{$chr}{$pos});
    print "$_\tna\t$clonal\n" unless(exists $hash{$chr}{$pos});
}
