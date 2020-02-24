#! bin/env -perl
$file=$ARGV[0]; # grep out SBW
open FILE, "$file" || die "no input file!";


while (<FILE>) {
	chomp $_;
	if ($_ =~ /##/) {next;}
	if ($_ =~ /#CHROM/) {@names=split(/\s+/,$_);
		foreach (@names[0,1,3,9..$#names]) {
			print $_." ";
		}
		print "\n";
		next;	}
	@line=split(/\s+/,$_);
	#print $line[0]." ".$line[1]." ".$line[3]." ".$line[4]." ";
	if ($line[4] =~ ",") {next;}
	if ($line[6] eq "my_snp_filter") {next;}
#		print $line[0]." ".$line[1]." ".$line[3]." ".$line[4]." ";
#	if ($line[7] =~ /MQ=(\d{2}\.\d{2})/) {$mq_snp=$1;}
#	if ($line[5]<100 or $mq_snp <59) { next;}
	print $line[0]." ".$line[1]." ".$line[3]." ";
	foreach $n (@line[9..$#line]) {
#	print $names[$counter]." ";
	@sample=split(/:/,$n);
		        if (($n =~ ":") and ($n !~ ",")) {
                               print $line[3]." ";
				}
                       if (($n =~ ":") and ($n =~ ",")) {
				@depth=split(/,/,$sample[1]);
				if (($depth[0]+$depth[1])== 0) {
                                         print "- ";
						next;	}
				if ($depth[0]/($depth[0]+$depth[1])>=0.99) {
      					  print $line[3]." ";}
      				if ($depth[1]/($depth[0]+$depth[1]) >= 0.99){
      					  print $line[4]." ";}
      				if (($depth[0]/($depth[0]+$depth[1])<0.99) and ($depth[0]/($depth[0]+$depth[1]) > 0.01) and ($depth[1]/($depth[0]+$depth[1]) < 0.99) and ($depth[1]/($depth[0]+$depth[1]) >0.01 )) {print "-"." ";}
				}
 #       		if ($n =~ /\./) {print "- ";}
	}
print "\n";
}
