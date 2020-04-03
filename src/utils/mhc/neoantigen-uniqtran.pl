#!/usr/bin/env perl 
use strict;
use warnings;
use utf8;

my $file = shift;
open IN, "<$file" or die "$!";
my %hash;
my %NMhash;
while(<IN>){
    chomp;
    if($_ =~ m/^>/){
        my($gene, $NM) = (split " ", $_)[1, 2];
        my $seq = <IN>;
        chomp($seq);
        next if($NM !~ m/NM/);
        $hash{$gene}{$NM} = $seq;
        $NMhash{$NM} = $gene;
    }
}
close IN;

$file = shift;
open IN, "$file" or die $!;
my %ID = ();
while(<IN>){
    chomp;
    next if($_ =~ /^hgncid/);
    my ($gencode, $ucsc) = (split /\t/, $_)[6, 8];
    $ID{$gencode} = $ucsc;
}
close IN;

my $file2 = shift;
open INF, "<$file2" or die "$!";
while(<INF>){
    chomp;
    if($_ =~ m/^PeptidePostion/){
        print "$_\tFlagTranscript\tTotalTranscript\n";
        next;
    }
    next if($_ =~ /missing/);
    my ($widetype, $gene, $NM) = (split "\t", $_)[7, 8, 11];
    if($NM !~ /^NM/){
        $NM =~ s/\_.+$//;
        if(exists $ID{$NM}){
            $NM = $ID{$NM};
        }
    }
    unless(exists $hash{$gene}{$NM}){
        if(exists $NMhash{$NM}){
            $gene = $NMhash{$NM};
        }
    }
    my $num_NM = 1;
    if(exists $hash{$gene}){
        $num_NM = scalar keys %{$hash{$gene}};
    }
    my $NM_n = 0;
    foreach my $NM_i(keys %{$hash{$gene}}){
        if($hash{$gene}{$NM_i} =~ m/$widetype/){
            $NM_n ++;
        }
    }
    $NM_n = 1 if($NM_n == 0);
    print "$_\t$NM_n\t$num_NM\n";
}

#===============================================================================
#                         This is a lovely split line ~                        
#===============================================================================                         
#                                                                              
#         FILE: selectuniqtranscript.pl                                                     
#                                                                              
#        USAGE: ./selectuniqtranscript.pl                                                   
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
#      CREATED: 2017年08月15日 14时47分26秒                                                  
#     REVISION: ---                                                            
#==============================================================================


