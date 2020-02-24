#!/usr/bin/perl

my $posFile = $ARGV[0];
my $joinList = $ARGV[1];
my $outFile = $ARGV[2];

open (POS, "$posFile") || die "cant load file";
open (LIST, "$joinList") || die "cant load file";
open (DEST, ">$outFile") || die "cant load file";

$| = 1;
%posHash = ();
while(<POS>)
{
	chomp;
	@input = split(/\s+/, $_);
	$poskey = join(" ", $input[0], $input[1], $input[2]);
	$posHash{$poskey} = $poskey;
}

while(<LIST>)
{
	chomp;
	$fileName = $_;
	%hash = ();
	open (SUBS, "$fileName") || next;
	while (<SUBS>){
		chomp;
		@input = split(/\s+/, $_);
		$newkey = join(" ", $input[0], $input[1], $input[2]);
		$hash{$newkey} = join(" ", $input[4]+ $input[5]+ $input[6]+ $input[7]+ $input[8]+ $input[9]+ $input[10]+ $input[11]);#sum the value.join(" ", $input[4], $input[5], $input[6], $input[7], $input[8], $input[9], $input[10], $input[11])
	}
	foreach $key (keys %posHash){
		if ($hash{$key} eq undef){
			$posHash{$key} = join(" ", $posHash{$key},0,0,0,0,0,0,0,0);#if we dont have the position,we get zero.($posHash{$key},0,0,0,0,0,0,0,0)
		}
		else{
			$posHash{$key} = join(" ", $posHash{$key},$hash{$key});
		}
	}
	print $fileName, "\n";
}

@sorted = sort {(split /\s+/, $a)[0] cmp (split /\s+/, $b)[0] || (split /\s+/, $a)[1] <=> (split /\s+/, $b)[1]} keys %posHash;


foreach $element (@sorted){
	print DEST $posHash{$element}, "\n";
}

