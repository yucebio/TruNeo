if(@ARGV < 1)
{
	die "Usage: perl $0 '/*/*.hla.top,bwahla' '/*/*.hla.txt,polysolver' > hla1.result 2> hla2.result\n";
}

my %hash1 = ();
my %hash2 = ();
my $t1n = 0;
my $t2n = 0;
foreach my $i (0..$#ARGV)
{
	my ($file,$tool) = split /,/,$ARGV[$i];
	open IN,$file or die $!;
	if($tool =~ /bwahla/i)
	{
		$t1n++;
		$t2n++;
		while(<IN>)
		{
			chomp;
			my ($hla1,$hla2) = (split)[1,2];
			$hla1 =~ s/(HLA-\w+\*\d+:\d+):.*/$1/;
			$hla2 =~ s/(HLA-\w+\*\d+:\d+):.*/$1/;
			if($hla1 =~ /HLA-D/)
			{
				if(!exists $hash2{$hla1}{$tool})
				{
					$hash2{$hla1}{'t'}++;
					$hash2{$hla1}{$tool} = 1;
				}
				if(!exists $hash2{$hla2}{$tool})
				{
					$hash2{$hla2}{'t'}++;
					$hash2{$hla2}{$tool} = 1;
				}
			}
			else
			{
				$hla1 =~ s/\*//g;
				$hla2 =~ s/\*//g;
				if(!exists $hash1{$hla1}{$tool})
				{
					$hash1{$hla1}{'t'}++;
					$hash1{$hla1}{$tool} = 1;
				}
				if(!exists $hash1{$hla2}{$tool})
				{
					$hash1{$hla2}{'t'}++;
					$hash1{$hla2}{$tool} = 1;
				}
			}
		}
	}
	elsif($tool =~ /polysolver/i)
	{
		$t1n++;
		while(<IN>)
		{
			chomp;
			my @arr = split;
			my @tmp1 = split /_/,$arr[1];
			my @tmp2 = split /_/,$arr[2];
			if(!exists  $hash1{"$arr[0]$tmp1[2]:$tmp1[3]"}{$tool})
			{
				$hash1{"$arr[0]$tmp1[2]:$tmp1[3]"}{'t'}++;
				$hash1{"$arr[0]$tmp1[2]:$tmp1[3]"}{$tool} = 1;
			}
			if(!exists  $hash1{"$arr[0]$tmp2[2]:$tmp2[3]"}{$tool})
			{
				$hash1{"$arr[0]$tmp2[2]:$tmp2[3]"}{'t'}++;
				$hash1{"$arr[0]$tmp2[2]:$tmp2[3]"}{$tool} = 1;
			}
		}
	}
	close IN;
}

my @T1 = ();
my @T2 = ();
foreach my $i (sort keys %hash1)
{
	push @T1,$i if($hash1{$i}{'t'} == $t1n);
}
foreach my $i (sort keys %hash2)
{
	push @T2,$i if($hash2{$i}{'t'} == $t2n);
}

my $HLA1 = join ',', @T1;
print "$HLA1\n";
my $HLA2 = join ',', @T2;
print STDERR "$HLA2\n";


#===============================================================================
#
#         FILE: overlapHLA.pl
#
#        USAGE: ./overlapHLA.pl
#
#  DESCRIPTION:
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Longyun Chen, chenly@yucebio.com
# ORGANIZATION: YuceBio
#      VERSION: 1.0
#      CREATED: 2017年03月17日 10时57分21秒
#     REVISION: ---
#===============================================================================
