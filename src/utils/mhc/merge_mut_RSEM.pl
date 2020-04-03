#!/usr/bin/perl -w 
use strict;
use Getopt::Long;
my ($mut,$exp,$vaf,$out,$help);
GetOptions(
	"m=s"  => \$mut,
	"e=s"  => \$exp,
	"v=s"  => \$vaf,
	"o=s"  => \$out,
	"h"    => \$help,
);

my $usage=<<INFO;
Usage:
perl $0 [options]

Options:
	-m:The neoantigen annotation summary table.
	-e:The gene expression data of gencode.
	-v:The rna vaf.
	-o:the output file.
	-h:Get the usage.
INFO

die $usage unless($mut && $exp && $vaf && $out);
die $usage if ($help);

open MUT,$mut or die $!;
open RSEMEXP,$exp or die $!;
open VAF,$vaf or die $!;
open OUT,">$out" or die $!;
print OUT join("\t",("Chr","Pos","Ref","Alt","Gene","Effect","RNAdep","RNAvaf","RC","TPM","FPKM",)),"\n";

my %exp;
while(<RSEMEXP>)
{
	chomp;
	#gene_name       gene_id transcript_id(s)        length  effective_length        expected_count  TPM     FPKM
	my($rc,$tpm,$fpkm,$gene) = (split /\t/,$_)[4,5,6,7];
	$exp{$gene}{"tpm"}  = $tpm;
        $exp{$gene}{"fpkm"} = $fpkm;
       	$exp{$gene}{"rc"}   = $rc;
}
close RSEMEXP;

my %vafh = ();
while(<VAF>)
{
	chomp;
	my ($chr,$pos,$fre,$alt) = (split /\t/,$_)[0,1,2,4];
	$chr =~ s/^chr//;
	my $key = join "_", $chr, $pos;
	$vafh{$key}{"vaf"} = $fre;
	$vafh{$key}{"rnamut"} = $alt;
}
close VAF;

while(<MUT>)
{
	chomp;
	next if($_ =~ /^Chr/);
	my($chr,$pos,$ref,$alt,$gene,$effect) = (split /\t/,$_)[0,1,2,3,5,6];
	$chr =~ s/^chr//;
	my $flag = join("_",($chr,$pos));
	my $rnadep = 0;
	my $rnavaf = 0;
	if(exists $vafh{$flag})
	{
		$rnavaf = $vafh{$flag}{"vaf"};
		$rnadep = $vafh{$flag}{"rnamut"};
	}
	if(!exists $exp{$gene}{"tpm"})
	{
		$exp{$gene}{"tpm"} = 0;
	}
	if(!exists $exp{$gene}{"fpkm"})
	{
		$exp{$gene}{"fpkm"} = 0;
	}	
	if(!exists $exp{$gene}{"rc"})
	{
		$exp{$gene}{"rc"} = 0
	}
	my ($fpkm,$rc,$tpm) = (0) x 3;
	$fpkm = $exp{$gene}{"fpkm"};
	$tpm  = $exp{$gene}{"tpm"};
	$rc   = $exp{$gene}{"rc"};
	print OUT join("\t",($chr,$pos,$ref,$alt,$gene,$effect,$rnadep,$rnavaf,$rc,$tpm,$fpkm)),"\n";
}
close MUT;
