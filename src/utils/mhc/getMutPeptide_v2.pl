#!/usr/bin/env perl 
use strict;
use warnings;
use utf8;

our %codon1 = (TTT=>"F", TTC=>"F", TCT=>"S", TCC=>"S", TAT=>"Y", TAC=>"Y", TGT=>"C", TGC=>"C", TTA=>"L", TCA=>"S", TAA=>"*", TGA=>"*", TTG=>"L", TCG=>"S", TAG=>"*", TGG=>"W", CTT=>"L", CTC=>"L", CCT=>"P", CCC=>"P", CAT=>"H", CAC=>"H", CGT=>"R", CGC=>"R", CTA=>"L", CTG=>"L", CCA=>"P", CCG=>"P", CAA=>"Q", CAG=>"Q", CGA=>"R", CGG=>"R", ATT=>"I", ATC=>"I", ACT=>"T", ACC=>"T", AAT=>"N", AAC=>"N", AGT=>"S", AGC=>"S", ATA=>"I", ACA=>"T", AAA=>"K", AGA=>"R", ATG=>"M", ACG=>"T", AAG=>"K", AGG=>"R", GTT=>"V", GTC=>"V", GCT=>"A", GCC=>"A", GAT=>"D", GAC=>"D", GGT=>"G", GGC=>"G", GTA=>"V", GTG=>"V", GCA=>"A", GCG=>"A", GAA=>"E", GAG=>"E", GGA=>"G", GGG=>"G");

sub shortaa
{
    my %short = (
        "Gly"=>"G",
        "Ala"=>"A",
        "Val"=>"V",
        "Leu"=>"L",
        "Ile"=>"I",
        "Phe"=>"F",
        "Pro"=>"P",
        "Ser"=>"S",
        "Thr"=>"T",
        "His"=>"H",
        "Trp"=>"W",
        "Cys"=>"C",
        "Asp"=>"D",
        "Glu"=>"E",
        "Lys"=>"K",
        "Tyr"=>"Y",
        "Met"=>"M",
        "Asn"=>"N",
        "Gln"=>"Q",
        "Arg"=>"R"
    );
    my $rawaa = shift;
    my @aa = split //,$rawaa;
    my $newaa = '';
    for(my $i = 0; $i < @aa; $i++) {
        if($aa[$i] =~ /[A-Z]/) {
            my $tmp = "$aa[$i]$aa[$i+1]$aa[$i+2]";
            my $tmpn = $tmp;
            $tmpn = $short{$tmp} if(exists $short{$tmp});
            $newaa .= $tmpn;
            $i = $i + 2;
        } else {
            $newaa .= $aa[$i];
        }
    }
    return $newaa;
}

my %complementary = (
    "A"=>"T",
    "T"=>"A",
    "G"=>"C",
    "C"=>"G",
    "N"=>"N"
);

my ($file, $genefile, $fastafile, $dnafasta, $out, $germf, $phasfile) = @ARGV;

my %germsnp = ();
if (-e $germf) {
    open GERM, $germf or die $!;
    while (<GERM>) {
        chomp;
        my @arr = split /\t/, $_;
        next if($arr[0] eq "chrom");
        $arr[10] =~ s/%//;
        next if($arr[10] < 20);
        $germsnp{$arr[0]}{$arr[1]} = $arr[3];
    }
    close GERM;
}

my (%queue, %need_trans);
open GENE, "<$genefile" or die "$!";
my (%mrnastart, %mrnaend, %mrnastrand, %mrnasnp);
while (<GENE>) {
    s/[\r\n]+$//;
    next if($_ =~ /exonCount/);
    my @field = split (/\t/, $_);
    @field >= 11 or die "Error: invalid record found in gene file (>=11 fields expected): <$_>\n";
    shift @field;         #refGene and ensGene has bin as the first column

    my ($name, $chr, $strand, $txstart, $txend, $cdsstart, $cdsend, $exonstart, $exonend) = @field[0, 1, 2, 3, 4, 5, 6, 8, 9];
    #$need_trans{$name} or next;
    my ($mrnastart, $mrnaend);

    #next we need to make sure that there is no intron between transcription start and translation start (this is rare but it happens when cdsstart is not in the first exon)
    my @exonstart = split (/,/, $exonstart);
    my @exonend = split (/,/, $exonend);
    $txstart++;
    $cdsstart++;
    @exonstart = map {$_+1} @exonstart;

    $mrnastrand{$name} = $strand;
    if ($strand eq '+') {
        #<---->-----<--->----<------>----<----->-----<--->
        #             **********
        my $intron = 0;
        for my $i (0 .. @exonstart-1) {
            $i and $intron += ($exonstart[$i]-$exonend[$i-1]-1);
            if ($cdsstart >= $exonstart[$i] and $cdsstart <= $exonend[$i]) {
                $mrnastart = $cdsstart-$txstart+1-$intron;
            }
            if ($cdsend >= $exonstart[$i] and $cdsend <= $exonend[$i]) {
                $mrnaend = $cdsend-$txstart+1-$intron;
            }
            if (exists $germsnp{$chr}) {
                foreach my $pp (sort {$a<=>$b} keys %{$germsnp{$chr}}) {
                    last if($pp > $cdsend);
                    last if($pp > $exonend[$i]);
                    next if($pp < $cdsstart);
                    if ($i == 0 and $pp < $exonstart[$i]) {
                        delete $germsnp{$chr}{$pp};
                    }
                    if ($pp >= $exonstart[$i] and $pp <= $exonend[$i]) {
                        my $tranpp = $pp-$txstart+1-$intron;
                        $mrnasnp{$name}{$tranpp} = $germsnp{$chr}{$pp};
                    }
                }
            }
        }
    } elsif ($strand eq '-') {
        #<---->-----<--->----<------>----<----->-----<--->
        #             **********
        my $intron = 0;
        for (my $i=@exonstart-1; $i>=0; $i--) {
            $i<@exonstart-1 and $intron += ($exonstart[$i+1]-$exonend[$i]-1);
            if ($cdsend >= $exonstart[$i] and $cdsend <= $exonend[$i]) {
                $mrnastart = $txend-$cdsend+1-$intron;
            }
            if ($cdsstart >= $exonstart[$i] and $cdsstart <= $exonend[$i]) {
                $mrnaend = $txend-$cdsstart+1-$intron;
            }
            if (exists $germsnp{$chr}) {
                foreach my $pp (sort {$a<=>$b} keys %{$germsnp{$chr}}) {
                    last if($pp > $cdsend);
                    last if($pp > $exonend[$i]);
                    next if($pp < $cdsstart);
                    if ($i == 0 and $pp < $exonstart[$i]) {
                        delete $germsnp{$chr}{$pp};
                    }
                    if ($pp >= $exonstart[$i] and $pp <= $exonend[$i]) {
                        my $tranpp = $txend-$pp+1-$intron;
                        $mrnasnp{$name}{$tranpp} = $complementary{$germsnp{$chr}{$pp}};
                    }
                }
            }
        }
    }

    $mrnastart{$name} = $mrnastart;
    $mrnaend{$name} = $mrnaend;
}
close GENE;

my %phasing = ();
if (-e $phasfile) {
    open IN,$phasfile or die $!;
    while(<IN>) {
        chomp;
        my @arr = split /\t/,$_;
        my @pos = split /,/,$arr[2];
        next if(@pos < 2);
        my @alt = split /,/,$arr[3];
        $arr[1] =~ s/^chr//;
        foreach my $i (0..$#pos) {
            $phasing{"$arr[1]:$pos[$i]:$alt[$i]"} = $arr[0];
        }
    }
    close IN;
}

my %id = ();
my %cluster = ();
my %sort = ();
open IN, "<$file" or die "$!";
while(<IN>){
    chomp;
    my ($gene, $gc, $gs, $alt, $transcript, $type, $cchange, $Achange) = (split "\t", $_)[ 4, 0, 1, 3, 6, 5, 7, 8]; 
    $cchange =~ s/^c.//;
    $cchange =~ s/^n.//;
    #$transcript =~ s/\.\d+$//;
    $Achange = &shortaa($Achange);
    my ($start, $end, $obs, $spl1, $spl2) = (0, 0, 0, 0, 0);      
    if ($cchange =~ m/^(\d+)\w>(\w)$/) {
        ($start, $end, $obs) = ($1, $1, $2);
    } elsif ($cchange =~ m/^(\d+)_(\d+)delins(\w+)/) {
        #block substitution
        ($start, $end, $obs) = ($1, $2, $3);
    } elsif ($cchange =~ m/^(\d+)_(\d+)del\w+ins(\w+)/) {
        #block substitution
        ($start, $end, $obs) = ($1, $2, $3);
    } elsif ($cchange =~ m/^(\d+)del\w+ins(\w+)/) {
        #block substitution
        ($start, $end, $obs) = ($1, $1, $2);
    } elsif ($cchange =~ m/^(\d+)del\w+/) {  
        #single base deletion
        ($start, $end, $obs) = ($1, $1, '');  
    } elsif ($cchange =~ m/^(\d+)_(\d+)del(\w*)/) { 
        #multi-base deletion 
        ($start, $end, $obs) = ($1, $2, '');
    } elsif ($cchange =~ m/^(\d+)_(\d+)ins(\w+)/) { 
        #insertion
        ($start, $end, $obs) = ($2, $2, $3); 
        #print "$start, $end, $obs\n";
        #if end is equal to start, this is an insertion 
    } elsif ($cchange =~ m/^(\d+)_(\d+)dup(\w+)/) { 
        #insertion
        my $inspos = $2 + 1;
        ($start, $end, $obs) = ($inspos, $inspos, $3); 
    } elsif ($cchange =~ m/^(\d+)dup(\w+)/) { 
        #dup
        my $inspos = $1 + 1;
        ($start, $end, $obs) = ($inspos, $inspos, $2);
    } elsif ($cchange =~ m/^(\d+)_(\d+)(\w+)/) {
        #non-frameshift substitution 
        ($start, $end, $obs) = ($1, $2, $3); 
    } else { 
        if($type =~ m/spli/){
            next;
            if($cchange =~ m/\d+-\d+/){
                next;
            }else{
                if($cchange =~ m/^(\d+)\+(\d+)\w>(\w)/){
                    ($start, $end, $obs, $spl1, $spl2) = ($1, $1, $3, $2, $2);
                }elsif($cchange =~ m/^(\d+)\+(\d+)del\w+/){
                    ($start, $end, $obs, $spl1, $spl2) = ($1, $1, '', $2, $2);
                }elsif($cchange =~ m/^(\d+)\+(\d+)_(\d+)\+(\d+)del(\w+)/){
                    ($start, $end, $spl1, $spl2, $obs) = ($1, $3, $2, $4, '');
                }elsif($cchange =~ m/^(\d+)\+(\d+)_(\d+)\+(\d+)ins(\w+)/){
                    ($start, $end, $spl1, $spl2, $obs) = ($1, $3, $2, $4, $5);
                }elsif($cchange =~ m/^(\d+)\+(\d+)ins(\w+)/){
                    ($start, $end, $obs, $spl1, $spl2) = ($1, $1, $3, $2, $2);
                }else{
                    next;
                }
                $gs = $gs - $spl1 if($mrnastrand{$transcript} eq "+");
                $gs = $gs + $spl1 if($mrnastrand{$transcript} eq "-");
            }
        }
        #   die "Error: invalid coding change format: <$cchange> within <$_>\n";
    }
    my $key = "$gene:$transcript";
    if(!exists $id{$key}) {
        $id{$key} = 0;
    }
    if(!exists $cluster{$key}){
        $cluster{$key} = 0;
    }
    if(!exists $queue{"$key:$id{$key}"} || !exists $phasing{"$gc:$gs:$alt"} || $phasing{"$gc:$gs:$alt"} != $cluster{$key}) {
        $id{$key}++;
        $cluster{$key} = $phasing{"$gc:$gs:$alt"} if(exists $phasing{"$gc:$gs:$alt"});
    }
    push @{$queue{"$key:$id{$key}"}}, [$start, $end, $spl1, $spl2, $obs, $cchange, $type, $gc, $gs, $Achange];
    $sort{"$key:$id{$key}"}{'chr'} = $gc;
    $sort{"$key:$id{$key}"}{'pos'} = $gs;
    $need_trans{$transcript}++;
}
close IN;

open OUTD, ">$out.details" or die "$!";
open OUTR, ">$out.result.txt" or die "$!";
open OUTW, ">$out.widetype.txt" or die "$!";
open OUTH, ">$out.peptideID.txt" or die "$!";
open OUTE, ">$out.conflict.txt" or die "$!";
print OUTE "Gene\tIsoform\tAnnodb\tWildaa\tMutaa\n";

open FASTA, "<$fastafile" or die "$!";
my (%mrnaseq);
my ($curname, $curseq);
while (<FASTA>) {
    s/[\r\n]+$//; 
    if (m/^>([\w\.]+)/) {
        if ($curseq) {
            $mrnaseq{$curname} = $curseq;
        }
        $curname = $1;
        $curseq = '';
    } else { 
        $curseq .= $_;
    }
    $curseq and $mrnaseq{$curname} = $curseq; #process the last sequence
}
close FASTA;

open FASTAD, "<$dnafasta" or die "$!";
my(%dnaseq);
my $curentchr = "";
my $seq = "";
while(<FASTAD>){
    chomp($_);
    if ($_ =~ m/>/) {
        if($curentchr){
            $dnaseq{$curentchr} = $seq;
        }
        $curentchr = $_;
        $curentchr = $1 if(($curentchr =~ m/>chr([X|Y])/ || $curentchr =~ m/>chr(\d+)/) && $curentchr !~ m/_/);
        print "$curentchr#######\n";
        $seq = '';
    }else{
        $seq .= $_;
    }
}
close FASTAD;
$dnaseq{$curentchr} = $seq;
foreach my $chr(keys %dnaseq){
    print "#################$chr\n";
}

my $peptideID = 0;
foreach my $i (sort {$sort{$a}{'chr'} cmp $sort{$b}{'chr'} or $sort{$a}{'pos'}<=>$sort{$b}{'pos'} or $a cmp $b} keys %queue){
    my ($line, $transcript, $id) = split /:/,$i;
    print "$line, $transcript, $id\n";
    #print "$line, $transcript, $start, $end, $obs, $cchange, $type, , $gc, $gs,$spl1, $spl2\n";
    #$verbose and print STDERR "NOTICE: Processing $line with $transcript, $start, $end, $obs, $cchange\n";
    if (not defined $mrnaseq{$transcript}) {
        print STDERR "WARNING: cannot find mRNA sequence for $transcript in the fastafile $fastafile\n";
        next;
    }
    if (not defined $mrnastart{$transcript}) {  
        print STDERR "WARNING: cannot find annotation for $transcript in the genefile $genefile or cannot infer the transcription start site\n";
        next; 
    }
    if (not defined $mrnaend{$transcript}) {
        print STDERR "WARNING: cannot find annotation for $transcript in the genefile $genefile or cannot infer the transcription end site\n";  
        next; 
    }
    my $dna = substr ($mrnaseq{$transcript}, $mrnastart{$transcript}-1, $mrnaend{$transcript}-$mrnastart{$transcript}+1) ;
    my @dna = split (//, $dna);               
    my ($protein1, $protein2);
    my $warning = '';
    if (exists $mrnasnp{$transcript}) {
        foreach my $pp (keys %{$mrnasnp{$transcript}}) {
            $dna[$pp - $mrnastart{$transcript}] = $mrnasnp{$transcript}{$pp};
        }
        $dna = join ('', @dna);
    }
    $protein1 = translateDNA ($dna); 
    next unless($protein1 ne "");
    print "$dna\n$protein1\n";
    if($protein1 =~ m/\*\w/){
        #i don't know why, but it works~
        $dna = substr ($mrnaseq{$transcript}, $mrnastart{$transcript}, $mrnaend{$transcript}-$mrnastart{$transcript}+2);
        @dna = split (//, $dna);
        if (exists $mrnasnp{$transcript}) {
            foreach my $pp (keys %{$mrnasnp{$transcript}}) {
                $dna[$pp - $mrnastart{$transcript} - 1] = $mrnasnp{$transcript}{$pp};
            }
            $dna = join ('', @dna);
        }
        $protein1 = translateDNA ($dna);
        if($protein1 =~ m/\*\w/){
            $dna = substr ($mrnaseq{$transcript}, $mrnastart{$transcript}+1, $mrnaend{$transcript}-$mrnastart{$transcript}+3) ;
            @dna = split (//, $dna);
            if (exists $mrnasnp{$transcript}) {
                foreach my $pp (keys %{$mrnasnp{$transcript}}) {
                    $dna[$pp - $mrnastart{$transcript} - 2] = $mrnasnp{$transcript}{$pp};
                }
                $dna = join ('', @dna);
            }
            $protein1 = translateDNA ($dna);
            if($protein1 =~ m/\*\w/){
                $dna = substr ($mrnaseq{$transcript}, $mrnastart{$transcript}-1, $mrnaend{$transcript}-$mrnastart{$transcript}+1);
                @dna = split (//, $dna);
                if (exists $mrnasnp{$transcript}) {
                    foreach my $pp (keys %{$mrnasnp{$transcript}}) {
                        $dna[$pp - $mrnastart{$transcript}] = $mrnasnp{$transcript}{$pp};
                    }
                    $dna = join ('', @dna);
                }
                $protein1 = translateDNA ($dna);
            }
        }
    }
    print "$mrnastart{$transcript},$mrnaend{$transcript},$dna\n";
    my $wtdna = $dna;
    my ($function, $aachange, $cDNA, $pchange, $indelobs) = ('') x 5;
    my ($gchr, $gpos) = ('', '');
    my $stopflag = 0;
    my @aapos = ();
    foreach my $j (0..$#{$queue{$i}}) {
        my ($start, $end, $spl1, $spl2, $obs, $cchange, $type, $gc, $gs, $Achange) = @{$queue{$i}[$j]};
        if ($end > length ($mrnaseq{$transcript})) {
            print STDERR "ERROR: transcript end ($mrnaend{$transcript}) for $transcript is longer than transcript length ${\(length ($mrnaseq{$transcript}))}, skipping this transcript\n";
            $protein2 = "";
            last;
        } else {
            if ($end > @dna) {
                print STDERR "ERROR in $line: end position of variant ($end) in $transcript is longer than coding portion length ${\(scalar @dna)}, skipping this transcript\n";
                $protein2 = "";
                last;
            }
            if(($spl1 && $spl2) || $type =~ m/stop_lost/){
                my $pre = substr($dna, 0, $start);
                my $late = "";
                if($mrnastrand{$transcript} eq "+"){
                    $late = substr($dnaseq{$gc}, $gs, 90);
                }else{
                    $late = substr($dnaseq{$gc}, $gs-90, 90);
                    $late = reverse($late);
                    $late =~ tr/actg/tgac/;
                }
                print "dna $dna\n";
                $dna = $pre.$late;
                $dna = lc($dna);
                $start = $start + $spl1;
                $end = $end + $spl2;
                #print "$dna\n";
            }
            @dna = split (//, $dna); 
            if ($start == $end and $cchange !~ m/ins/ and $cchange !~ m/dup/ and $cchange !~ m/del/) {
                #this is a substition
                $dna[$start - 1] = $obs;
            } elsif ($cchange =~ m/\d+del\w+ins\w+/) {
                #this is a multiple substition
                splice (@dna, $start-1, $end-$start+1, $obs);
                #my @obstmp = split //,$obs;
                #foreach my $i (0..$#obstmp) {
                #    my $p = $start + $i - 1;
                #    $dna[$p] = $obstmp[$i];
                #}
            } elsif ($start == $end and $obs) {  
                #this is an insertion
                splice (@dna, $start-1, 0, $obs);
            } else { 
                #this is a deletion
                splice (@dna, $start-1, $end-$start+1, $obs);
            }
            $dna = join ('', @dna);
            #print "$dna\n";
            my $aastart = int(($start-1)/3)+1; 
            my $aaend1 = int (($end-1)/3)+1;
            $aaend1 = int (($start+length($obs)-1)/3)+1 if(($cchange =~ m/ins/ or $cchange =~ m/dup/) and $cchange !~ /\d+del\w+ins\w+/);
            my $aaend2 = $aastart;
            if($cchange =~ /\d+del\w+ins\w+/) {
                $aaend2 = int (($start+length($obs)-1-1)/3)+1;
            } elsif ($cchange =~ m/ins/ or $cchange =~ m/dup/) {
                $aaend2 = int (($start+length($obs)-1)/3)+1;
            }
            $protein2 = translateDNA ($dna);
            if ($stopflag > 0 && $aastart >= $stopflag) {
                $protein2 = "";
            }
            #print "$protein1\n$protein2\n";
            last unless($protein2 ne "");
            if (substr ($protein2, $aastart-1, $aaend2-$aastart+1) =~ m/\*/) {
                $stopflag = $aastart;
                if (length($function) == 0) {
                    $function = 'immediate-stopgain';
                } else {
                    $function .= ';immediate-stopgain';
                }
            } elsif ($end>=$mrnaend{$transcript}-2) {
                if (length($function) == 0) {
                    $function = 'immediate-stoploss';
                } else {
                    $function .= ';immediate-stoploss';
                }
            } else {
                if ($protein1 eq $protein2) {
                    $function = 'synonymous';
                } else {
                    if (length($function) == 0) {
                        $function = $type;
                    } else {
                        $function .= ";$type";
                    }
                }
            }
            my $wtaa = substr ($protein1, $aastart-1, $aaend1-$aastart+1);
            my $mtaa = substr ($protein2, $aastart-1, $aaend2-$aastart+1);
            my @wt = split //,$protein1;
            my @mt = split //,$protein2;
            my $index = $aastart-1;
            my $wtchange = 'na';
            my $change = 'na';
            if($Achange =~ /^p\.([A-Z]+)(\d+)(.+)$/) {
                $wtchange = $1;
                $index = $2-1;
                $change = $3;
            }
            if($Achange =~ /ins|del|dup|fs/) {
                if($wt[$index] ne $wtchange) {
                    print OUTE "$line\t$transcript\t$Achange\t$wt[$index]\t-\n";
                    $protein2 = "";
                    last;
                }
            } elsif($Achange =~ /\*$/) {
                if($wt[$index] ne $wtchange || $mt[$index] ne '*') {
                    print OUTE "$line\t$transcript\t$Achange\t$wt[$index]\t$mt[$index]\n";
                    $protein2 = "";
                    last;
                }
            } elsif($Achange =~ /\?$/) {
                if($wt[$index] ne $wtchange) {
                    print OUTE "$line\t$transcript\t$Achange\t$wt[$index]\t$mt[$index]\n";
                    $protein2 = "";
                    last;
                }
            } elsif($wt[$index] ne $wtchange || $mt[$index] ne $change) {
                print OUTE "$line\t$transcript\t$Achange\t$wt[$index]\t$mt[$index]\n";
                $protein2 = "";
                last;
            }
            if (length($aachange) == 0) {
                $aachange = "(position $aastart-$aaend1 changed from " . $wtaa . ' to ' . $mtaa;
            } else {
                $aachange .= ";position $aastart-$aaend1 changed from " . $wtaa . ' to ' . $mtaa;
            }
            if (length($cDNA) == 0) {
                $cDNA = "c.$cchange";
            } else {
                $cDNA .= ";c.$cchange";
            }
            if (length($pchange) == 0) {
                $pchange = $Achange;
            } else {
                $pchange .= ";$Achange";
            }
            $gchr = $gc if(length($gchr) == 0);
            if (length($gpos) == 0) {
                $gpos = $gs;
            } else {
                $gpos .= ";$gs";
            }
            push @aapos, $aastart;
            $indelobs = $obs;
        }
    }
    next unless($protein2 ne "");
    $aachange .= ')';
    $protein2 =~ s/\*.+/\*/;
    #delete anything after the STOP codon 
    my $subaa = "";
    my $subaawide = "";
    my $subdna = "";
    my $subdnawt = "";
    my $pepstart = (sort {$a<=>$b} @aapos)[0];
    my $pepend = (sort {$b<=>$a} @aapos)[0];
    my $s = ($pepstart-15 < 0)?0:$pepstart-15;
    my $l = $pepstart-$s;
    my $wll2 = ($l+$pepend-$pepstart+14+$s > length($protein1))? length($protein1) - $s:$l+$pepend-$pepstart+14;
    $subaawide = substr($protein1, $s, $wll2);
    $subdnawt = substr($wtdna, $s*3, $wll2*3);
    if($function =~ m/stopgain/){
        $subaa = substr($protein2, $s, $l+$pepend-$pepstart);
        $subdna = substr($dna, $s*3, ($l+$pepend-$pepstart)*3);
    }else{
        my $ll2 = ($l+$pepend-$pepstart+14+$s > length($protein2))? length($protein2) - $s:$l+$pepend-$pepstart+14;
        if($cDNA =~ m/ins|dup/ && length($indelobs)%3!=0 && $cDNA !~ m/\d+_\d+del\w+ins\w+/){
            $ll2 = length($protein2) - $s;
        } elsif ($cDNA =~ m/\d+del([ATGC]+)$/) {
            if (length($1)%3 != 0) {
                $ll2 = length($protein2) - $s;
            }
        }
        print "$transcript, $protein1, $protein2, $s, $l, $ll2\n";
        if($subaa = substr($protein2, $s, $ll2)){
            $subdna = substr($dna, $s*3, $ll2*3);
        }else{
            next;
        }
    }
    $subdnawt =~ s/(.{100})/$1\n/g;
    $subdna =~ s/(.{100})/$1\n/g;
    $wtdna =~ s/(.{100})/$1\n/g;
    $dna =~ s/(.{100})/$1\n/g;
    $protein2 =~ s/(.{100})/$1\n/g; 
    $protein1 =~ s/(.{100})/$1\n/g;
    $peptideID ++;
    $pchange = "" unless($pchange);
    my $rank = '';
    foreach my $p (@aapos) {
        my $r = $p-$pepstart+$l;
        if (length($rank) == 0) {
            $rank = $r;
        } else {
            $rank .= ";$r";
        }
    }
    my $peprideinfo = join "\t", $peptideID, $line, $gchr, $gpos, $transcript, $cDNA, $function, $aachange, $pchange, $subaa, $rank;
    print OUTH "$peprideinfo\n";
    print OUTD ">Peptide$peptideID $line $transcript $cDNA WTsubcdna\n";
    print OUTD "$subdnawt\n";
    print OUTD ">Peptide$peptideID $line $transcript $cDNA MTsubcdna\n";
    print OUTD "$subdna\n";
    print OUTD ">Peptide$peptideID $line $transcript $cDNA WTsubpeptide\n";
    print OUTD "$subaawide\n";
    print OUTD ">Peptide$peptideID $line $transcript $cDNA MTsubpeptide\n";
    print OUTD "$subaa\n";
    print OUTD ">Peptide$peptideID $line $transcript $cDNA WTcdna\n";
    print OUTD "$wtdna\n";
    print OUTD ">Peptide$peptideID $line $transcript $cDNA MTcdna\n";
    print OUTD "$dna\n";
    print OUTD ">Peptide$peptideID $line $transcript $cDNA WTpeptide\n";
    print OUTD "$protein1\n";
    print OUTD ">Peptide$peptideID $line $transcript $cDNA MTpeptide\n";
    print OUTD "$protein2\n";
    next if($function eq 'synonymous' || $function eq 'immediate-stopgain' || $warning =~ m/error/);
    print OUTW ">Peptide$peptideID\n";
    print OUTW "$subaawide\n";
    print OUTR ">Peptide$peptideID\n";
    print OUTR "$subaa\n";
}

sub translateDNA{
    my ($seq) = @_; 
    my ($nt3, $protein); 
    $seq = uc $seq;  
    #length ($seq) % 3 == 0 or printerr "WARNING: length of DNA sequence to be translated is not multiples of 3: <length=${\(length $seq)}>\n"; 
    while ($seq =~ m/(...)/g) {       
        #defined $codon1{$1} or print "WARNING: invalid triplets found in DNA sequence to be translated: <$1> in <$seq>\n" and die;
        defined $codon1{$1} or return "";
        $protein .= $codon1{$1}; 
    }
    return $protein;
}

#===============================================================================
                                                                                                                                                                                                                                                                                                                                                                                                                                                            #
#         FILE: getpeptide.pl
                                                                                                                                                                                                                                                                                                                                                                                                                                                            #
#        USAGE: ./getpeptide.pl  
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
#      CREATED: 2015年09月29日 15时59分13秒
#     REVISION: ---
#===============================================================================


