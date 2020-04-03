#!/usr/bin/env perl
use strict;
use warnings;
use utf8;

my ($file,$outd) = @ARGV;
open IN, "<$file" or die "$!";
while(<IN>){
    chomp;
    if($_ =~ m/>/){
        my $peptide = $_;
        $peptide =~ s/>//;
        `echo "$_" > $outd/$peptide.fa`;
        my $seq = <IN>;
        chomp($seq);
        `echo $seq >> $outd/$peptide.fa`;
        print "$outd/$peptide.fa\n";
    }
}

#===============================================================================
#                         This is a lovely split line ~                        
#===============================================================================                         
#                                                                              
#         FILE: generatenetctlpanformat.pl                                                     
#                                                                              
#        USAGE: ./generatenetctlpanformat.pl                                                   
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
#      CREATED: 2017年08月29日 15时37分05秒                                                  
#     REVISION: ---                                                            
#==============================================================================


