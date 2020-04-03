#!/usr/bin/env perl 
use strict;
use warnings;
use utf8;

my ($fa, $file, $i) = @ARGV;
die $0 if(@ARGV < 3);
#my $fa = "/mnt/cfs/prj15a/M25C000/pep_alignment/database/Homo_sapiens.GRCh37.75.pep.all.fa";
#my $fa = "/mnt/nfs/user/wangjq/MyDatabase/database/protein/Homo_sapiens.GRCh38.pep.all.fa.gz";
if($fa =~ /\.gz$/)
{
	open FASTA,"gzip -dc $fa |" or die $!;
}
else
{
	open FASTA,$fa or die "$!";
}
my (%mrnaseq);
my ($curname, $curseq);
while (<FASTA>) {
    chomp;
    if (m/^>(.*)/) {
        if ($curseq && length($curseq) > 7) {
            $mrnaseq{$curname} = $curseq;
        }
        $curname = $1; 
        $curseq = '';
    } else {
        $curseq .= $_;   
    }   
}
close FASTA;
if ($curseq && length($curseq) > 7) {
    $mrnaseq{$curname} = $curseq; #process the last sequence
}

open IN, "<$file" or die "$!";
while(<IN>){
    chomp;
    if($_ =~ m/Peptide/){
        print "$_\tSelfilter\n";
        next;
    }
    my $peptide = (split "\t", $_)[$i];
    my $flag = 0;
    foreach my $id(keys %mrnaseq){
        if($mrnaseq{$id} =~ m/$peptide/){
            $flag = 1;
            last;
        }
    }
    print "$_\t$flag\n";
}
close IN;

#===============================================================================
#
#         FILE: uniq.pl
#
#        USAGE: ./uniq.pl  
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
#      CREATED: 2016年05月24日 21时27分58秒
#     REVISION: ---
#===============================================================================


