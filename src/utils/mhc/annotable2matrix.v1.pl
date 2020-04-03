#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use Data::Dumper;
use utf8;
use Encode;

my ($annotable,$type,$outable1,$outable2,$help);
GetOptions(
	"a=s"  => \$annotable,
	"t=s"  => \$type,
	"o1=s" => \$outable1,
	"o2=s" => \$outable2,
    "h"    => \$help,
);

my $usage=<<INFO;
Usage:
perl $0 [options]

Options:
	-a:The neoantigen annotation summary table.
	-t:The type of sequency.DNA or RNA.
	-o1:The neoantigen MHC I/II type table.
	-o2:The High scoe neoantigen MHC I/II type table.
	-h:Get the usage.
INFO

die $usage unless($annotable);
die $usage if ($help);

$type ||= "DNA";

my @title1;
my @title2;
if($type eq "DNA")
{
	@title1=("NeoRank","Score","HLA","Peptide","NeoepitopeScore","NeoepitopeType","WildtypeScore","WildtypePeptide","Gene","Transcript","CDSMutation","Function","ProteinChangeAbbr","ExtendedPeptide","TotalTranscript","Selfilter","AnchorFlag","WildPeptideMHCPrediction","WildTAPPredictionScore","WildCleavagePredictionScore","WildCombinedPredictionScore","PeptideMHCPrediction","TAPPredictionScore","CleavagePredictionScore","CombinedPredictionScore","CosmicCount","KnownNeoepitopeDBFlag","DNAVAF","CCLEexpStatus","CCF","Clonal","PeptideLength");
	@title2=("NeoRank","Score","HLA","Peptide","NeoepitopeScore","NeoepitopeType","WildtypeScore","WildtypePeptide","Gene","Transcript","CDSMutation","Function","ProteinChangeAbbr","ExtendedPeptide","TotalTranscript","Selfilter","CosmicCount","DNAVAF","CCLEexpStatus","CCF","Clonal","PeptideLength");
}
else 
{
	@title1=("NeoRank","Score","HLA","Peptide","NeoepitopeScore","NeoepitopeType","WildtypeScore","WildtypePeptide","Gene","Transcript","CDSMutation","Function","ProteinChangeAbbr","ExtendedPeptide","TotalTranscript","Selfilter","AnchorFlag","WildPeptideMHCPrediction","WildTAPPredictionScore","WildCleavagePredictionScore","WildCombinedPredictionScore","PeptideMHCPrediction","TAPPredictionScore","CleavagePredictionScore","CombinedPredictionScore","CosmicCount","KnownNeoepitopeDBFlag","DNAVAF","ExpStatus","RNAdepth","RNAVAF_TPM","CCF","Clonal","PeptideLength");
	@title2=("NeoRank","Score","HLA","Peptide","NeoepitopeScore","NeoepitopeType","WildtypeScore","WildtypePeptide","Gene","Transcript","CDSMutation","Function","ProteinChangeAbbr","ExtendedPeptide","TotalTranscript","Selfilter","CosmicCount","DNAVAF","ExpStatus","RNAdepth","RNAVAF_TPM","CCF","Clonal","PeptideLength");
}


my @otitle;
if($annotable =~ /\.anchor.TAP\./)
{
	@otitle = @title1;
}
else
{
	@otitle = @title2;
}
open OUT1,">$outable1" or die $!;
open OUT2,">$outable2" or die $!;

print OUT1 join("\t",@otitle),"\n";
print OUT2 join("\t",@otitle),"\n";

my %index;
open ANNOTABLE,$annotable or die $!;
while(<ANNOTABLE>)
{
	chomp;
	my @oline;
	if($_ =~ /^NeoRank/)
	{
		my @title = split /\t/,$_;
		for(my $i=0;$i<@title;$i++)
		{
			$index{$title[$i]} = $i;
		}
	}
	else
	{
		my @line = (split /\t/,$_);
		my $length;
		foreach my $title(@otitle)
		{
			if(exists $index{$title})
			{
				my $id = $index{$title};
				if($title eq "Peptide")
				{
					$length = length($line[$id]);
				}
				if(($title eq "CCLEexpStatus") || ($title eq "ExpStatus"))
				{
					my $exp = $line[$id];
					if($exp == 0 )
					{
						$line[$id] = "none";
					}
					elsif($exp == 1 )
					{
						$line[$id] = "low";
					}
					elsif($exp == 2)
					{
						$line[$id] = "medium";
					}
					elsif($exp == 3)
					{
						$line[$id] = "high";
					}
					else
					{
						print STDERR "Need check line:$_\n";
					}
				}
				push @oline,$line[$id];
			}
			else
			{
				next;
			}
		}
		push @oline,$length;
		print OUT1 join("\t",@oline),"\n";
		my $NeoepitopeType = $oline[5];
		if($NeoepitopeType =~/SB|WB/)
		{
			print OUT2 join("\t",@oline),"\n";
		}
	}
}
close ANNOTABLE;
close OUT1;
close OUT2;
