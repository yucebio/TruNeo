#!/usr/bin/env perl 
use strict;
use warnings;
use utf8;

my $file = shift;
open IN, "<$file" or die "$!";
my( $startflag, $endflag, $Affinity_i) = (0, 0, 0);
my %hash;
while(<IN>){
    chomp;
    if($_ =~ m/Affinity\(nM\)/){
        $startflag = 1;
        $_ =~ s/\s+/\t/g;
        my @temp = split "\t", $_;
        foreach my $i(0..@temp-1){
            next unless($temp[$i] =~ m/Affinity\(nM\)/);
            $Affinity_i = $i;
            last;
        }
        next;
    }
    if($startflag && $_ =~ m/-----------/){
        $endflag ++;
        if($endflag == 2){
            $startflag = 0;
            $endflag = 0;
        }
        next;
    }
    next if($startflag == 0);
    $_ =~ s/\s+/\t/g;
    my($pos, $HLA, $peptide, $Identity, $Affinity) = (split "\t", $_)[1, 2, 3, 4, $Affinity_i];
    my $gene = (split "_", $Identity)[0];
    my $length = length($peptide);
    my $key = join ",", $pos, $HLA, $gene, $length;
    $hash{$key}{'aff'} = $Affinity;
    $hash{$key}{'pep'} = $peptide;
}

my $pepfile = shift;
my %position = ();
my %fram = ();
open IN, $pepfile or die $!;
while(<IN>){
    chomp;
    my @arr = split /\t/,$_;
    my @mult = split /;/,$arr[10];
    foreach my $i (@mult){
        push @{$position{"Peptide$arr[0]"}}, $i;
    }
    if($arr[6] =~ /frameshift/){
        $fram{"Peptide$arr[0]"} = 1;
    }
}
close IN;

my $Tfile = shift;
open INT, "<$Tfile" or die "$!";
while(<INT>){
    chomp;
    if($_ =~ m/Affinity\(nM\)/){
        $startflag = 1;
        next;
    }
    if($startflag && $_ =~ m/-----------/){
        $endflag ++;
        if($endflag == 2){
            $startflag = 0;
            $endflag = 0;
        }
        next;
    }
    next if($startflag == 0);
    $_ =~ s/\s+/\t/g;
    my($pos, $HLA, $peptide, $Identity, $Affinity) = (split "\t", $_)[1, 2, 3, 4, $Affinity_i];
    my $type = "";
    $type = "WB" if($Affinity < 500);
    $type = "SB" if($Affinity < 50);
    my $length = length($peptide);
    my $judge = 0;
    if(exists $position{$Identity}){
        foreach my $i (@{$position{$Identity}}){
            if($length + $pos >= $i && $pos < $i){
                $judge = 1;
                last;
            }
        }
    } else{
        $judge = 1 if($length + $pos >= 15 && $pos < 15);
    }
    next unless ($judge);
    my $gene = (split "_", $Identity)[0];
    my $key = join ",", $pos, $HLA, $gene, $length;
    my $widetype = "";
    $widetype = $hash{$key}{'aff'} if(exists $hash{$key});
    my $widetypepep = "";
    $widetypepep = $hash{$key}{'pep'} if(exists $hash{$key});
    if(exists $fram{$Identity} && $widetype eq ""){
        my $postmp = $pos;
        while($postmp > 0){
            $postmp--;
            $key = join ",", $postmp, $HLA, $gene, $length;
            if(exists $hash{$key}){
                $widetype = $hash{$key}{'aff'};
                $widetypepep = $hash{$key}{'pep'};
                last;
            }
        }
    }
    if($widetype eq "" ){
        next;
    }
    my $outpos = $pos + 1;
    my $out = join "\t", $outpos, $HLA, $peptide, $Identity, $Affinity, $type, $widetype, $widetypepep;
    print "$out\n";
}
#===============================================================================
#
#         FILE: CAPneoantigen-08-netMHCpan_select.pl
#
#        USAGE: ./CAPneoantigen-08-netMHCpan_select.pl  
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
#      CREATED: 2015年12月30日 11时49分53秒
#     REVISION: ---
#===============================================================================


