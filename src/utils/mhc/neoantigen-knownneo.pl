#!/usr/bin/env perl 
use strict;
use warnings;
use utf8;

#Level 1: validated successfully with matched HLA
#Level 2: validated successfully with no matched HLA
#Level 3: select but no validation with matched HLA
#Level 4: select but no validation with no matched hla
#Level 5: validated failure

my $file = shift;
open IN, "<$file" or die "$!";
my %hash;
my %PeptideHash;
while(<IN>){
    chomp;
    next if($_ =~ m/^#/);
    my($hla, $Peptide, $Epitope, $Info) = (split "\t", $_)[2, 3, 4, 6];
    my $Level = 0;
    $hla =~ s/([A|B|C])(\d+:\d+)/HLA-$1*$2/ if($hla && $hla ne "");
    if($Info =~ m/=1/ && $hla && $hla ne "" && $hla =~ m/:/){
        $Level = "1";
        $hash{$Epitope}{$hla} = $Level;
        $PeptideHash{$Peptide} = $Level;
    }elsif($Info =~ m/=1/){
        $Level = "2";
        $hash{$Epitope}{'na'} = $Level;
        $PeptideHash{$Peptide} = $Level;
    }elsif($Info =~ m/APC=0/ || $Info =~ m/PRR=0;POR=0/){
        $Level = 5;
        $hash{$Epitope}{'na'} = $Level;
    }
}

my $neoantigen = shift;
open INN, "<$neoantigen" or die "$!";
while(<INN>){
    chomp;
    if($_ =~ m/PeptidePostion/){
        print "$_\tKnownNeoepitopeDBFlag\n";
        next;
    }
    my($hla, $Epitope) = (split "\t")[1, 2];
    my $Level = 0;
    if(exists $hash{$Epitope}{$hla}){
        $Level = "$hash{$Epitope}{$hla}A";
    }elsif(exists $hash{$Epitope}{'na'}){
        $Level = "$hash{$Epitope}{'na'}A";
    }else{
        foreach my $Peptide(keys %PeptideHash){
            if($Peptide =~ m/$Epitope/){
                $Level = "$PeptideHash{$Peptide}B";
            }
        }
    }
    print "$_\t$Level\n";
}

#===============================================================================
#                         This is a lovely split line ~                        
#===============================================================================                         
#                                                                              
#         FILE: CAPneoantigen-11-KnownNeoepitope.pl                                                     
#                                                                              
#        USAGE: ./CAPneoantigen-11-KnownNeoepitope.pl                                                   
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
#      CREATED: 2017年08月30日 10时31分14秒                                                  
#     REVISION: ---                                                            
#==============================================================================


