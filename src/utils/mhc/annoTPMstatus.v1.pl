#!/usr/bin/env perl 
use strict;
use warnings;
use Data::Dumper;
use utf8;

my $in1=$ARGV[0];
### $in1 include rnavaf depth tpm file
my $in2=$ARGV[1];
### $in2 sample.transcript.fp.anchor.TAP.hotspot.KnownNeoepitope.dnavaf

my %hash;
open IN1, "$in1" or die "$!";
while(<IN1>){
    chomp;
    next if($_ =~ m/^Chr/);
    my($chr,$pos,$ref,$alt,$gene,$dep,$vaf,$tpm) = (split /\t/,$_)[0,1,2,3,4,6,7,9];
    next if($vaf == 0);
    $chr =~ s/^chr//;
    my $key = join "\t", $gene, $pos;
    $hash{$key}{"tpm"} = $vaf*$tpm;
    $hash{$key}{"depth"} = $dep;
}
close IN1;

my @rank;
my $i = scalar keys %hash;
foreach my $key(keys %hash){
    push @rank, $hash{$key}{"tpm"};
}
@rank = sort{$a<=>$b} @rank;

my $lower = $rank[int($i*0.25)];
my $upper = $rank[int($i*0.75)];

open IN2, $in2 or die "$!";
while(<IN2>){
    chomp;
    if($_ =~ m/^PeptidePostion/){
        print "$_\tExpStatus\tRNAdepth\tRNAVAF_TPM\n";
        next;
    }
    my $status = 0;
    my ($gene, $mut) = (split "\t", $_)[8, 10];
    my $key = join "\t", $gene, $mut;
    my $rnavaf = 0;
    my $tpm    = 0;
    if(!exists $hash{$key})
    {
        print "$_\t0\t0\t0\n";
    }
    else
    {
        my $alt_tpm = $hash{$key}{"tpm"};
        my $rnaread = $hash{$key}{"depth"};
        if($alt_tpm >= $upper){$status = 3;}
        elsif($alt_tpm == 0){$status = 0;}
        elsif($alt_tpm <= $lower){$status = 1;}
        else{$status = 2;}
        print "$_\t$status\t$rnaread\t$alt_tpm\n";
    }
}
close IN2;


#===============================================================================
#                         This is a lovely split line ~                        
#===============================================================================                         
#                                                                              
#         FILE: annoTPMstatus.pl                                                     
#                                                                              
#        USAGE: ./annoTPMstatus.pl                                                   
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
#      CREATED: 2017年08月15日 21时30分17秒                                                  
#     REVISION: ---                                                            
#==============================================================================


