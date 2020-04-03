#!/usr/bin/perl -w

use strict;
use File::Basename;
my($widetypelist, $mutlist)=@ARGV;
open WList,"$widetypelist";

my %hash;
while(my $sample=<WList>)
{
	chomp $sample;
	my $samplename=basename $sample;
	$samplename =~ s/Wild_//;
	$samplename =~ s/_fa//;
	my($peptideID, $hla, $length) = (split '\.', $samplename)[0, 1, 2];
	open Each,"$sample" or die "$!";
	while(my $lane=<Each>)
	{
		chomp $lane;
		next if($lane !~ m/^\d/);
		my $pos = (split "\t", $lane)[0];
		my @temp = split "\t", $lane;
		shift @temp;
		my $o = join "\t", @temp;
		$hash{$peptideID}{$hla}{$length}{$pos} = $o;
	}
	close Each;
}
close WList;

print "PeptideID\tHLA\tLength\tPeptidePos\tPeptide\tPeptideMHCPrediction\tTAPPredictionScore\tCleavagePredictionScore\tCombinedPredictionScore\tNetCTLpanRank\tWildPeptide\tWildPeptideMHCPrediction\tWildTAPPredictionScore\tWildCleavagePredictionScore\tWildCombinedPredictionScore\tWildNetCTLpanRank\n";

open MList,"<$mutlist" or die "$!";
while(my $sample=<MList>)
{
	chomp $sample;
	my $samplename=basename $sample;
	$samplename =~ s/Mut_//;
	$samplename =~ s/_fa//;
	my($peptideID, $hla, $length) = (split '\.', $samplename)[0, 1, 2];
	open Each,"$sample" or die "$!";
	while(my $lane=<Each>)
	{
		chomp $lane;
		next if($lane !~ m/^\d/);
		my $pos = (split "\t", $lane)[0];
		my $widetype = "na\tna\tna\tna\tna\tna\tna";
		$widetype = $hash{$peptideID}{$hla}{$length}{$pos} if(exists $hash{$peptideID}{$hla}{$length}{$pos});
		print "$peptideID\t$hla\t$length\t$lane\t$widetype\n";
	}
	close Each;
}
close MList;
