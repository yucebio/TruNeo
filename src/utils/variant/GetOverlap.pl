#!/usr/bin/perl -w
use strict;
my ($file1, $file2, $i) = @ARGV;
open IN,'<',"$file1" or die "$!";
my %hash;
while(<IN>){
    chomp;
    my $Factor = (split "\t",$_)[0];
    $hash{$Factor} = 1 unless((split "\t",$_)[1]);
    $hash{$Factor} = (split "\t",$_)[1] if((split "\t",$_)[1]);
}
close IN;

open IN,'<',"$file2" or die "$!";
while(<IN>){
    chomp;
    my $Factor = (split "\t",$_)[$i];
    print "$hash{$Factor}\t$_\n" if (exists $hash{$Factor}); 
    print "0\t$_\n" unless (exists $hash{$Factor}); 
}
