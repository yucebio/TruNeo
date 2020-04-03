#!/usr/bin/env perl 
use strict;
use warnings;
use utf8;
use Encode;
use File::Basename;
use Getopt::Long;
use Data::Dumper;
use utf8;
use Encode;

my ($in1,$in2,$type,$out1,$out2,$help);
GetOptions(
    "i1=s"  => \$in1,
    "i2=s"  => \$in2,
    "t=s"   => \$type,
    "o1=s"  => \$out1,
    "o2=s"  => \$out2,
    "h"     => \$help,
);

my $usage=<<INFO;
Usage:
perl $0 [options]

Options:
    -i1:The neoantigen MHC I/II type table.
    -i2:The High scoe neoantigen MHC I/II type table.
    -t:The type of sequency.DNA or RNA.
    -o1:The neoantigen HLA summary table .
    -o2:The neoantigen HLA summary WBSB table.
    -h:Get the usage.
INFO

die $usage unless($in1 && $in2);
die $usage if ($help);

$type ||= "DNA";

my %highqc;
open IN2, "<$in2" or die "$!";
while(<IN2>)
{
	chomp;
	my($mutpep,$gene,$aachange) = (split /\t/,$_)[3,8,12];
	my $flag = join("\t",($gene,$aachange,$mutpep));
	$highqc{$flag} = 1;
}
close IN2;

open IN1, "<$in1" or die "$!";
my %hash;
my %hash_fac2;
my %info ;

while(<IN1>){
    chomp;
	next if($_ =~ /^NeoRank/);
    my @temp = split "\t", $_;
    my $factor1 = join "\t", $temp[8], $temp[12], $temp[3];
    my $factor2 = $temp[2];
    $hash{$factor1}{$factor2}{'mut'} = $temp[4];
    $hash{$factor1}{$factor2}{'wild'} = $temp[6];
    $hash_fac2{$factor2} = 1 ;
    if($type eq "DNA")
    {
        $info{$factor1} = join"\t",$temp[-4], $temp[-3], $temp[-2],$temp[-1];
    }
    else
    {
        $info{$factor1} = join"\t",$temp[-6], $temp[-5], $temp[-4], $temp[-3], $temp[-2], $temp[-1];
    }
}
close IN1;

open OUT1,">$out1" or die $!;
open OUT2,">$out2" or die $!;
print OUT1 "HLA Type\t\t";
print OUT2 "HLA Type\t\t";
my @factor;
foreach my $factor2(sort keys %hash_fac2){
    print OUT1 "\t$factor2\t$factor2";
    print OUT2 "\t$factor2\t$factor2";
    push @factor,$factor2;
}
print OUT1 "\n";
print OUT2 "\n";
my $length = @factor;
my @title = ("Mutated","Wildtype") x $length;
my $title = join("\t",("Gene","Mutation","Neopeptide",@title));
#my $rnasupport = encode_utf8("RNA中变异reads支持数");
#my $rnasupport = "Mut_Reads_In_RNA";
#my $clonalstat = encode_utf8("克隆状态");
#my $peptidelen = encode_utf8("肽段长度");
if($type eq "DNA")
{
    print OUT1 join("\t",($title,"CCLEexpStatus","CCF","Clonal","PeptideLength")),"\n";
    print OUT2 join("\t",($title,"CCLEexpStatus","CCF","Clonal","PeptideLength")),"\n";
}
else
{
    print OUT1 join("\t",($title,"ExpStatus","RNAdepth","RNAVAF_TPM","CCF","Clonal","PeptideLength")),"\n";
    print OUT2 join("\t",($title,"ExpStatus","RNAdepth","RNAVAF_TPM","CCF","Clonal","PeptideLength")),"\n";
}

foreach my $factor1(sort keys %hash){
    print OUT1 "$factor1";
    foreach my $factor2(@factor){
		if(exists $hash{$factor1}{$factor2})
		{
			print OUT1 "\t$hash{$factor1}{$factor2}{'mut'}\t$hash{$factor1}{$factor2}{'wild'}";
		}
    }
    print OUT1 "\t$info{$factor1}\n";
	if(exists $highqc{$factor1})
	{
		print OUT2 "$factor1";
		foreach my $factor2(@factor){
			if(exists $hash{$factor1}{$factor2})
			{
				print OUT2 "\t$hash{$factor1}{$factor2}{'mut'}\t$hash{$factor1}{$factor2}{'wild'}";
			}
		}
		print OUT2 "\t$info{$factor1}\n";
	}
}
close OUT1;
close OUT2;
