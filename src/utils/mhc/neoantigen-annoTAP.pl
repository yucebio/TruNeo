#!/usr/bin/perl -w
use strict;
my ($file1, $file2) = @ARGV;
open IN,'<',"$file1" or die "$!";
my %hash;
while(<IN>){
    chomp;
    my($peptideID, $HLA, $peptide, $Wildtype_MHC_Prediction, $Wildtype_TAP_Prediction_Score, $Wildtype_Cleavage_Prediction_Score, $Wildtype_Combined_Prediction_Score, $Mutation_MHC_Prediction, $Mutation_TAP_Prediction_Score,   $Mutation_Cleavage_Prediction_Score, $Mutation_Combined_Prediction_Score) = (split "\t", $_)[0, 1, 4, 11, 12, 13, 14, 5, 6, 7, 8];
    my $Factor = join "\t", $HLA, $peptide, $peptideID;
    $hash{$Factor} = join "\t", $Wildtype_MHC_Prediction, $Wildtype_TAP_Prediction_Score, $Wildtype_Cleavage_Prediction_Score, $Wildtype_Combined_Prediction_Score, $Mutation_MHC_Prediction, $Mutation_TAP_Prediction_Score, $Mutation_Cleavage_Prediction_Score, $Mutation_Combined_Prediction_Score;
}
close IN;

open IN,'<',"$file2" or die "$!";
while(<IN>){
    chomp;
    my($HLA, $peptide, $peptideID) = (split "\t", $_)[1, 2, 3];
    $HLA =~ s/\*//;
    my $Factor = join "\t", $HLA, $peptide, $peptideID;
    if(exists $hash{$Factor}){
        print "$_\t$hash{$Factor}\n";
    }else{
        print "$_\t0\t0\t0\t0\t0\t0\t0\t0\n";
    }
}
