#!/usr/bin/perl

my $snpFile = $ARGV[0];

open (SNP, "$snpFile") || die "cant load file";


$| = 1;
while (<SNP>)
{
	chomp;
	@input = split(/\s+/, $_);	
	$scaffold = $input[0];
	$position = $input[1];
	$ref = $input[2];
	$cons = $input[scalar(@input)-1];
	$lines=0;
	$diff=0;
	$change=0;
	for ($x=3;$x<(scalar(@input)-1);$x++){
		if ($input[$x] =~ /-/){
			next;
		}
		if ($input[$x] ne $cons){
			$change = $x;
			$diff++;	
		}
		$lines++;
	}
	if (($lines > 1) and ($diff < $lines) and ($diff > 0)){
		print $change-3, " ", $_, "\n";
	}
	
}



