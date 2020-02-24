#!/usr/bin/perl

my $snpFile = $ARGV[0];
my $consensusFile = $ARGV[1];

open (SNP, "$snpFile") || die "cant load file";
open (DEST, ">$consensusFile") || die "cant load file";

$| = 1;
while (<SNP>)
{
	chomp;
	@input = split(" ", $_);	
	$scaffold = $input[0];
	$position = $input[1];
	$ref = $input[2];

	print DEST join(" ",$scaffold,$position,$ref);

	$total = 0;
	$totalA = 0;
	$totalC = 0;
	$totalG = 0;
	$totalT = 0;
	#countNucs
	for ($x=3;$x<scalar(@input);$x+=8){
		$A=$input[$x];
		$a=$input[$x+1];
		$C=$input[$x+2];
		$c=$input[$x+3];
		$G=$input[$x+4];
		$g=$input[$x+5];
		$T=$input[$x+6];
		$t=$input[$x+7];

		$lineTotal=$A+$a+$C+$c+$G+$g+$T+$t;
		$totalA+=$A+$a;
		$totalC+=$C+$c;
		$totalG+=$G+$g;
		$totalT+=$T+$t;
		$total+=$lineTotal;

		if ($lineTotal eq 0){
			print DEST " -";
			next;
		}
		if ((($A+$a)/$lineTotal)>.80){
			if (($A>2) and ($a>2))
			{
				print DEST " A";
				next;
			}
		}
                if ((($C+$c)/$lineTotal)>.80){
                        if (($C>2) and ($c>2))
                        {
                                print DEST " C";
                                next;
                        }
                }
		if ((($G+$g)/$lineTotal)>.80){
                        if (($G>2) and ($g>2))
                        {
                                print DEST " G";
                                next;
                        }
                }
		if ((($T+$t)/$lineTotal)>.80){
                        if (($T>2) and ($t>2))
                        {
                                print DEST " T";
                                next;
                        }
                }
		print DEST " -";
	}
	if ($total eq 0){
		print DEST " -\n";
        	next;
        }
	if (($totalA/$total)>.50){
		print DEST " A";
	}
        if (($totalC/$total)>.50){
                print DEST " C";
        }
        if (($totalG/$total)>.50){
                print DEST " G";
        }
        if (($totalT/$total)>.50){
                print DEST " T";
        }
	print DEST "\n";
}

