#!/usr/bin/env perl 
use strict;
use warnings;
use utf8;

my $anno = shift;
open IN, "<$anno" or die "$!";
my %trans;
while(<IN>){
    chomp;
    next if($_ =~ m/^Chr/);
    my @arr = split "\t", $_;
    my $tr = $arr[7];
    #if($arr[7] =~ /^([^\.]+)\.\d+$/){
    #    $tr = $1;
    #}
    my $key = join "\t", $arr[5], $arr[0], $arr[1], $tr;
    $trans{$key} = 1;
}

my $file = shift;
open IN, "<$file" or die "$!";
my %sort;
while(<IN>){
    chomp;
    if($_ =~ m/^PeptidePostion/){
        print "NeoRank\tClassfication\tScore\tScore_m\t$_\n";
        next;
    }
    my @arr = split "\t", $_;
    my($pos, $HLA, $peptide, $peptideID, $affinity, $type, $widescore, $widepeptide, $gene, $chr, $pos_genome, $NM, $cdsmut, $muttype, $AAchange, $PPchange, $peptide_raw, $changepos, $flag_transcript, $total_transcript, $selfilter, $anchor_flag, $Psite, $Wildtype_MHC_Prediction, $Wildtype_TAP_Prediction_Score, $Wildtype_Cleavage_Prediction_Score, $Wildtype_Combined_Prediction_Score, $Mutation_MHC_Prediction, $Mutation_TAP_Prediction_Score, $Mutation_Cleavage_Prediction_Score, $Mutation_Combined_Prediction_Score, $COSMIC_COUNT, $KnownNeoepitopeDBFlag, $vaf, $CCLE_status) = split "\t", $_;
    next if($muttype eq "stop_lost");
    $CCLE_status = $arr[-1];
    my $key = join "\t", $pos, $HLA, $peptide, $gene, $chr, $pos_genome, $CCLE_status;
    my $pos_genome1 = (split /;/,$pos_genome)[0];
    my $key2 = join "\t", $gene, $chr, $pos_genome1, $NM;
    #my $sigma = 1 - log($affinity)/log(50000);
    my $sigma = &tanh((500 - $affinity) / 200);
    my $sigma2 = &tanh(($widescore / $affinity - 1) / 2) / 2 + 0.5;
    if($total_transcript == 0){
        $total_transcript = 1;
        $flag_transcript = 1;
    }
    $Mutation_TAP_Prediction_Score = &tanh($Mutation_TAP_Prediction_Score * 2.5) / 2 + 0.5;
    $Mutation_Cleavage_Prediction_Score = &tanh($Mutation_Cleavage_Prediction_Score * 3);
    my $combine = 0.8 * $sigma + 0.15 * $Mutation_Cleavage_Prediction_Score + 0.05 * $Mutation_TAP_Prediction_Score;
    $combine = 0 if($combine < 0);
    #if($CCLE_status == 0){
    #    $CCLE_status = 1;
    #}
    if(@arr > 35){
        $CCLE_status = 1 - 0.7 ** $CCLE_status;
    }
    if($CCLE_status == 0){
        $CCLE_status = 0.2;
    }
    my $score = $combine * $sigma2 * $CCLE_status * (1 - $selfilter) * $flag_transcript * $vaf/$total_transcript;
    $score = 0 if($score < 0);
    my $normal = 27/100;
    if($type ne ""){
        if($muttype =~ m/^fram/ && $type eq "HP"){
            my $finalscore = $score * 16/$normal;
            if(!exists $sort{$key} || $sort{$key}{'tr'} == 0){
                $sort{$key}{'info'} = $_;
                $sort{$key}{'final'} = $finalscore;
                $sort{$key}{'class'} = 1;
                $sort{$key}{'score'} = $score;
            }
        }elsif($anchor_flag > 0 && $type eq "HP"){
            my $finalscore = $score * 9.6/$normal;
            if(!exists $sort{$key} || $sort{$key}{'tr'} == 0){
                $sort{$key}{'info'} = $_;
                $sort{$key}{'final'} = $finalscore;
                $sort{$key}{'class'} = 2;
                $sort{$key}{'score'} = $score;
            }
        }elsif($type eq "HP"){
            my $finalscore = $score * 8/$normal;
            if(!exists $sort{$key} || $sort{$key}{'tr'} == 0){
                $sort{$key}{'info'} = $_;
                $sort{$key}{'final'} = $finalscore;
                $sort{$key}{'class'} = 3;
                $sort{$key}{'score'} = $score;
            }
        }elsif($muttype =~ m/^fram/){
            my $finalscore = $score * 4/$normal;
            if(!exists $sort{$key} || $sort{$key}{'tr'} == 0){
                $sort{$key}{'info'} = $_;
                $sort{$key}{'final'} = $finalscore;
                $sort{$key}{'class'} = 4;
                $sort{$key}{'score'} = $score;
            }
        }elsif($anchor_flag > 0){
            my $finalscore = $score * 2.4/$normal;
            if(!exists $sort{$key} || $sort{$key}{'tr'} == 0){
                $sort{$key}{'info'} = $_;
                $sort{$key}{'final'} = $finalscore;
                $sort{$key}{'class'} = 5;
                $sort{$key}{'score'} = $score;
            }
        }else{
            my $finalscore = $score * 2/$normal;
            if(!exists $sort{$key} || $sort{$key}{'tr'} == 0){
                $sort{$key}{'info'} = $_;
                $sort{$key}{'final'} = $finalscore;
                $sort{$key}{'class'} = 6;
                $sort{$key}{'score'} = $score;
            }
        }
    }else{
        my $finalscore = $score * 0.2/$normal;
        if(!exists $sort{$key} || $sort{$key}{'tr'} == 0){
            $sort{$key}{'info'} = $_;
            $sort{$key}{'final'} = $finalscore;
            $sort{$key}{'class'} = 7;
            $sort{$key}{'score'} = $score;
        }
    }
    if(exists $trans{$key2}){
        $sort{$key}{'tr'} = 1;
    }elsif(!exists $sort{$key}{'tr'}){
        $sort{$key}{'tr'} = 0;
    }
}

my $rank = 0;
foreach my $key(sort {$sort{$b}{'final'}<=>$sort{$a}{'final'} or $sort{$a}{'class'}<=>$sort{$b}{'class'} or $sort{$b}{'score'}<=>$sort{$a}{'score'} or $a cmp $b} keys %sort){
    my $finalscore = $sort{$key}{'final'};
    my $i = $sort{$key}{'class'};
    my $score = $sort{$key}{'score'};
    $rank ++;
#   next if($finalscore == 0);
    my $finalscorep = $finalscore;
    print "$rank\t$i\t";
    printf "%.2f",$finalscorep;
    print "\t$score\t$sort{$key}{'info'}\n";     
}

sub tanh{
    my $value = shift;
    my $cal = 1;
    if($value < -100){
        $cal = -1;
    }elsif($value < 100){
        $cal = (exp($value) - exp(-$value)) / (exp($value) + exp(-$value));
    }
    return $cal;
}

#===============================================================================
#                         This is a lovely split line ~                        
#===============================================================================                         
#                                                                              
#         FILE: sort.pl                                                     
#                                                                              
#        USAGE: ./sort.pl                                                   
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
#      CREATED: 2017年08月15日 15时11分37秒                                                 
#     MODIFIED: 2019.6.5 by D.X, danxu@yucebio.com
#     REVISION: ---                                                            
#==============================================================================


