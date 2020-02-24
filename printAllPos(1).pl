#!/usr/bin/perl
my $infile = $ARGV[0];
open (SOURCE, "$infile") || die "cant load file";

$/ = ">";
$| = 1;	

while(<SOURCE>)
{
	chomp;
    	my ($header, $sequence) = split(/\n/,$_,2);
	next unless ($header && $sequence);
	my ($seqid) = $header =~ /^(\S+)/;
	$sequence =~ s/[\n\r]//g;   			
	@input = split("", $sequence);
	$n=1;	
	foreach $element (@input)
	{
		print $seqid, " ", $n, " ", uc($element), "\n";
		$n++;
	}
}

