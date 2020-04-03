#!/usr/bin/env perl 
use strict;
use warnings;
use utf8;

my ($hotspot, $file) = @ARGV;
#my $hotspot = "/mnt/cfs/prj17a/M373000/bin/hotspotbin/recurrent.position.stat.txt";
open IN, "<$hotspot" or die "$!";
my %hash;
while(<IN>){
    chomp;
    my($chr, $position, $count) = (split "\t", $_)[0, 1, 3];
    $hash{$chr}{$position} = $count;
}
close IN;

open INF, "<$file" or die "$!";
while(<INF>){
    chomp;
    if($_ =~ m/^PeptidePostion/){
        print "$_\tCosmicCount\n";
        next;
    }
    my($chr, $position) = (split "\t", $_)[9, 10];
    my $count = 0;
    my @pos = split /;/,$position;
    foreach my $i (@pos){
        $count = $hash{$chr}{$i} if(exists $hash{$chr}{$i});
    }
    print "$_\t$count\n";
}

#===============================================================================
#                         This is a lovely split line ~                        
#===============================================================================                         
#                                                                              
#         FILE: annohotspot.pl                                                     
#                                                                              
#        USAGE: ./annohotspot.pl                                                   
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
#      CREATED: 2017年08月16日 20时56分00秒                                                  
#     REVISION: ---                                                            
#==============================================================================


