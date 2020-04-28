#!/usr/bin/perl -w

open(PLAN,"<$ARGV[0]") || die;
my %hndl;

$chr="chr0";
$idx = 0;
while(<PLAN>){
	chomp;
	@plan = split(/\t/);
	$plan[2] = "X" if ($plan[2] == 23);
	if("chr$plan[2]" ne $chr){
		close(VDATA) if($chr ne "chr0");
		$chr = "chr$plan[2]";
		open(VDATA,"gzip -dc 1000g_slim_related/data.merged.$chr.txt.gz |") || die;
		$idx = 0;
	}
	while($idx<$plan[3]) {
		$idx++;
		$data = <VDATA>;
		chomp($data);
	}
	@temp=split(//,$data);
	if(!exists($hndl{"$plan[0]_$plan[1]"})){
		open($hndl{"$plan[0]_$plan[1]"},"|gzip > data_sample/sample.$plan[0].$plan[1].txt.gz") || die;
	}
	print {$hndl{"$plan[0]_$plan[1]"}} join("\t",@temp[0..$#temp])."\n";
}
close(VDATA);
close(PLAN);
foreach $plan (keys %hndl){
	close($hndl{$plan}); 
}
