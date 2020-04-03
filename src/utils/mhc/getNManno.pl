#!/usr/bin/env perl 
use strict;
use warnings;
use utf8;

my $file = shift;
open IN, "<$file" or die "$!";
while(<IN>){
    chomp;
    my($Chr, $Pos, $Ref, $Alt, $Impact, $Gene, $Effect, $Transcript, $Biotype, $cDNA, $pp, $Isoforms) = (split "\t", $_)[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11];
    #if(!($Effect eq "downstream_gene_variant" || $Effect eq "upstream_gene_variant" || $Effect eq "intergenic_region" || $Effect eq "intron_variant" || $Effect eq "synonymous_variant" || $Effect =~ m/UTR/ || $Effect eq "" || !($pp) || $pp eq "")){
        #print "$Chr\t$Pos\t$Ref\t$Alt\t$Gene\t$Effect\t$Transcript\t$cDNA\t$pp\n";
    #}
    if($Impact =~ /MODERATE|HIGH/){
        if($Isoforms =~ m/,/){
            my @temp = split ",", $Isoforms;
            foreach my $Iso(@temp){
                my($Effect_i, $Impact_i, $Gene_i, $Transcript_i, $Biotype_i, $cDNA_i, $pp_i) = (split "\\|", $Iso)[1, 2, 3, 6, 7, 9, 10];
                #if(!($Effect_i eq "downstream_gene_variant" || $Effect_i eq "upstream_gene_variant" || $Effect_i eq "intergenic_region" || $Effect_i eq "intron_variant" || $Effect_i eq "synonymous_variant" || $Effect_i =~ m/UTR/ || $Effect_i eq "" || $Effect_i eq "non_coding_exon_variant" || !($pp_i) || $pp_i eq "" || $Effect_i eq "splice_region_variant&intron_variant" || $Effect_i eq "splice_acceptor_variant&intron_variant" || $Effect_i eq 'splice_acceptor_variant&intron_variant' || $Effect_i =~ m/^spli/)){
                if($Impact_i =~ /MODERATE|HIGH/ && $pp_i && $Biotype_i =~ /coding/i){
                    if($Biotype_i eq "Coding"){
                        $Transcript_i =~ s/\..+//;
                    }
                    print "$Chr\t$Pos\t$Ref\t$Alt\t$Gene_i\t$Effect_i\t$Transcript_i\t$cDNA_i\t$pp_i\n";
                }
            }
        }
        elsif($pp && $Biotype =~ /coding/i){
            if($Biotype eq "Coding"){
                $Transcript =~ s/\..+//;
            }
            print "$Chr\t$Pos\t$Ref\t$Alt\t$Gene\t$Effect\t$Transcript\t$cDNA\t$pp\n";
        }
    }
}

#===============================================================================
#
#         FILE: getNManno.pl
#
#        USAGE: ./getNManno.pl  
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
#      CREATED: 2015年09月29日 15时35分26秒
#     REVISION: ---
#===============================================================================


