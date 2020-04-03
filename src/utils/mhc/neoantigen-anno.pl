#!/usr/bin/env perl 
use strict;
use warnings;
use utf8;

my $peptideIDfile = shift;
open IN, "<$peptideIDfile" or die "$!";
my %hash;
while(<IN>){
    chomp;
    my($id, @temp) = split "\t", $_;
    $hash{"Peptide$id"} = join "\t", @temp;
}
close IN;

my $MHCresultfile = shift;
open INR, "<$MHCresultfile" or die "$!";
my $head = join "\t", "PeptidePostion", "HLA", "Peptide", "PeptideID", "NeoepitopeScore", "NeoepitopeType", "WildtypeScore", "WildtypePeptide", "Gene", "Chr", "Position", "Transcript", "CDSMutation", "Function", "ProteinChange", "ProteinChangeAbbr", "ExtendedPeptide", "AAchangePosition";

print "$head\n";
while(<INR>){
    chomp;
    my $PeptideID = (split "\t", $_)[3];
    print "$_\t$hash{$PeptideID}\n" if(exists $hash{$PeptideID});
    print "$_\tmissing\n" unless(exists $hash{$PeptideID});
}


#===============================================================================
#
#         FILE: CAPneoantigen-08-anno.pl
#
#        USAGE: ./CAPneoantigen-08-anno.pl  
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
#      CREATED: 2016年06月15日 15时38分36秒
#     REVISION: ---
#===============================================================================


