my ($list,$dat,$outp)=@ARGV;

die "Usage: perl $0 HLA-DQA1*02:01,HLA-DQA1*03:02 pseudosequences.dat OUTPUT.prefix\n" if(@ARGV < 3);

my @HLA = split /,/,$list;
my %hash = ();
foreach my $i (@HLA)
{
	my ($type1, $type2, $type3) = ($1, $2, $3) if($i =~ m/(D.*)\*(\d+):(\d+)/);
	my $str = "$type2$type3";
	$hash{$type1}{$str} = 1;
}

foreach my $DRB (keys %{$hash{'DRB1'}})
{
	`cut -f 1 $dat | grep DRB1 | grep $DRB >> $outp.MHCII.list`;
}

foreach my $DPA (keys %{$hash{'DPA1'}})
{
	foreach my $DPB (keys %{$hash{'DPB1'}})
	{
		`cut -f 1 $dat | grep DPA1$DPA-DPB1$DPB >> $outp.MHCII.list`;
	}
}

foreach my $DQA (keys %{$hash{'DQA1'}})
{
	foreach my $DQB (keys %{$hash{'DQB1'}})
	{
		`cut -f 1 $dat | grep DQA1$DQA-DQB1$DQB >> $outp.MHCII.list`;
	}
}

`sort $outp.MHCII.list | uniq > $outp.MHCII.list.temp`;
`mv $outp.MHCII.list.temp $outp.MHCII.list`;
