#!/usr/bin/perl -w

$filepath=$ARGV[0];
$filepath=~/variants\.(chr.*?)\.txt\.gz/;
$chr = $1;
open(VAR1,"gzip -dc 1000g_slim/variants.$chr.txt.gz |") || die;
open(VAR2,"gzip -dc $filepath |") || die;
open(DATA1,"gzip -dc 1000g_slim/data.$chr.txt.gz |") || die;
open(DATA2,"gzip -dc 1000g_slim_related/data.$chr.txt.gz |") || die;
open(VDATA1,"|gzip > data.merged.$chr.txt.gz") || die;
open(VDATA2,"|gzip > data.relate.$chr.txt.gz") || die;

$var1=<VAR1>;
$var2=<VAR2>;
$data1=<DATA1>;
$data2=<DATA2>;
$var1_num = 0; $var2_num = 0; $var_same = 0;
while(1) {
	last if(!$var1);
	chomp($data1);
	@temp1 = split(/\t/,$var1);
	@temp2 = split(/\t/,$var2) if($var2);
	if($var2 && $temp1[1] eq $temp2[1] && $temp1[3] eq $temp2[3] && $temp1[4] eq $temp2[4]){
		print VDATA1 $data1.$data2;
		print VDATA2 $data2;
		$var_same++;
		$var1=<VAR1>;
		$var2=<VAR2>;
		$data1=<DATA1>;
		$data2=<DATA2>;
	}elsif(!$var2 || $temp1[1]<$temp2[1]){
		print VDATA1 $data1."0000000000000000000000000000000\n";
		print VDATA2 "0000000000000000000000000000000\n";
		$var1_num++;
		$var1=<VAR1>;
		$data1=<DATA1>;
	}else{
		$var2=<VAR2>;
		$data2=<DATA2>;
		$var2_num++;
	}
}
close(VAR1);
close(VAR2);
close(DATA1);
close(DATA2);
close(VDATA1);
close(VDATA2);
print "$var_same\t$var1_num\t$var2_num\n";
