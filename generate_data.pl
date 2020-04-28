#!/usr/bin/perl -w

$filepath=$ARGV[0];
$regionfile=$ARGV[1];
open(REGIONFILE,"<$regionfile") || die;
my %region;
while(<REGIONFILE>) {
	chomp;
	@temp=split(/\t/);
	push(@{$region{$temp[0]}},[$temp[1],$temp[2]]);
}
close(REGIONFILE);
foreach $chr ( keys %region){
	$region{$chr}=[sort{$a->[0] <=> $b->[0] || $a->[1] <=> $b->[1]} @{$region{$chr}}];
}
$filepath=~/ALL\.(chr.*?)\./;
$chr = $1;
open(VCF,"gzip -dc $filepath |") || die;
open(VINFO,"|gzip > variants.$chr.txt.gz") || die;
open(VDATA,"|gzip > data.$chr.txt.gz") || die;
$cds_idx=0;

while(<VCF>) {
	chomp;
	next if(/^#/);
	@temp=split(/\t/);
	@temp4 = split(/,/,$temp[4]);
	@temp7 = split(/;/,$temp[7]);
	$AF = 0; $EAS_AF = 0; $AMR_AF = 0; $AFR_AF = 0; $EUR_AF = 0; $SAS_AF = 0;
	$CDS = 0;
	while (exists($region{$chr}) && $cds_idx<scalar(@{$region{$chr}}) && $region{$chr}->[$cds_idx][1]<$temp[1]){
		$cds_idx++;
	}
	if(exists($region{$chr}) && $cds_idx<scalar(@{$region{$chr}}) && $temp[1]>$region{$chr}->[$cds_idx][0] && $temp[1] <= $region{$chr}->[$cds_idx][1]){
		$CDS = 1;
	}
	for ($i = 0 ; $i < @temp7 ; $i++){
		@af = split(/=/,$temp7[$i]);
		if($af[0] eq "AF"){
			$AF = $af[1];
		}elsif($af[0] eq "EAS_AF"){
			$EAS_AF = $af[1];
		}elsif($af[0] eq "AMR_AF"){
			$AMR_AF = $af[1];
		}elsif($af[0] eq "AFR_AF"){
			$AFR_AF = $af[1];
		}elsif($af[0] eq "EUR_AF"){
			$EUR_AF = $af[1];
		}elsif($af[0] eq "SAS_AF"){
			$SAS_AF = $af[1];
		}
	}
	@AF = split(/,/,$AF);
	@EAS_AF = split(/,/,$EAS_AF);
	@AMR_AF = split(/,/,$AMR_AF);
	@AFR_AF = split(/,/,$AFR_AF);
	@EUR_AF = split(/,/,$EUR_AF);
	@SAS_AF = split(/,/,$SAS_AF);

	for ($i = 0 ; $i < @temp4 ; $i++){
		if($AF=~/,/){
			print VINFO "$temp[0]\t$temp[1]\t$temp[2]\t$temp[3]\t$temp4[$i]\t$AF[$i]\t$EAS_AF[$i]\t$AMR_AF[$i]\t$AFR_AF[$i]\t$EUR_AF[$i]\t$SAS_AF[$i]\t$CDS\n";
		}else{
			print VINFO "$temp[0]\t$temp[1]\t$temp[2]\t$temp[3]\t$temp4[$i]\t$AF\t$EAS_AF\t$AMR_AF\t$AFR_AF\t$EUR_AF\t$SAS_AF\t$CDS\n";
		}
		$data = '';
		for($j = 9 ; $j < @temp ; $j++){
			@gt = split(/\||\//,$temp[$j]);
			$num = 0;
			foreach $gt (@gt){
				$gt = substr($gt,0,index($gt,':')) if($gt=~/:/);
				$num++ if($gt ne '.' && $gt == $i+1);
			}
			$data.=$num;
		}
		print VDATA $data."\n";
	}
}
close(VCF);
close(VDATA);
close(VINFO);
