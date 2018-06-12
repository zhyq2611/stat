#!/usr/bin/perl -w
use strict;
my ($fqstat,$duplexlog,$mpeliup,$outtype) = @ARGV;
open FL,"$fqstat" or die $!;
$/ = "\n\n";
my $read1_raw = <FL>;
my $read1_clean = <FL>;
my $read2_raw = <FL>;
my $read2_clean = <FL>;
$read1_raw =~ /total reads: (\d+)\ntotal bases: (\d+)\nQ20 bases: (\d+)\(.*\)\nQ30 bases: (\d+)\((.*)\)/;
my $read1_raw_read = $1;
my $read1_raw_base = $2;
my $read1_raw_Q20_base = $3;
my $read1_raw_Q30_base = $4;
my $Raw_R1_Q30 = $5;
$read2_raw =~ /total reads: (\d+)\ntotal bases: (\d+)\nQ20 bases: (\d+)\(.*\)\nQ30 bases: (\d+)\((.*)\)/;
my $read2_raw_read = $1;
my $read2_raw_base = $2;
my $read2_raw_Q20_base = $3;
my $read2_raw_Q30_base = $4;
my $Raw_R2_Q30 = $5;
$read1_clean =~ /total reads: (\d+)\ntotal bases: (\d+)\nQ20 bases: (\d+)\(.*\)\nQ30 bases: (\d+)\(.*\)/;
my $read1_clean_read = $1;
my $read1_clean_base = $2;
my $read1_clean_Q20_base = $3;
my $read1_clean_Q30_base = $4;
$read2_clean =~ /total reads: (\d+)\ntotal bases: (\d+)\nQ20 bases: (\d+)\(.*\)\nQ30 bases: (\d+)\(.*\)/;
my $read2_clean_read = $1;
my $read2_clean_base = $2;
my $read2_clean_Q20_base = $3;
my $read2_clean_Q30_base = $4;
close FL;
$/ = "\n";
my $raw_read = sprintf"%.2f",($read1_raw_read + $read2_raw_read)/1000000;
my $raw_base = sprintf"%.2f",($read1_raw_base + $read2_raw_base)/1000000;
my $clean_read = sprintf"%.2f",($read1_clean_read + $read2_clean_read)/1000000;
my $clean_base = sprintf"%.2f",($read1_clean_base + $read2_clean_base)/1000000;
my $clean_rate_read = sprintf"%.2f",$clean_read*100/$raw_read;
my $clean_rate_base = sprintf"%.2f",$clean_base*100/$raw_base;
print "Raw read: $raw_read(M)\n";
print "Raw base: $raw_base(M)\n";
print "Raw R1 Q30: $Raw_R1_Q30\n";
print "Raw R2 Q30: $Raw_R2_Q30\n";
print "Clean read: $clean_read(M)\n";
print "Clean base: $clean_base(M)\n";
print "Clean rate(read): $clean_rate_read(%)\n";
print "Clean rate(base): $clean_rate_base(%)\n";

my $map_read = 0;
my $map_molecule = 0;
my $molecule_with_dup = 0;

open FL,"$duplexlog" or die $!;
while(<FL>){
	chomp;
	if($_ =~ /Mapping read = (\d+)/){
		$map_read = $1;
	}elsif($_ =~ /Mapping molecule = (\d+)/){
		$map_molecule = $1;
	}elsif($_ =~ /molecule with PCRdup = (\d+)/){
		$molecule_with_dup += $1;
	}
}
my $map_rate = sprintf"%.2f",$map_read*100/($read1_clean_read + $read2_clean_read);
my $dup_rate = sprintf"%.2f",$map_molecule*2*100/$map_read;
$map_read = sprintf"%.2f",$map_read/1000000;
$molecule_with_dup = sprintf"%.2f",$molecule_with_dup/1000000;
$map_molecule = sprintf"%.2f",$map_molecule/1000000;
print "MAP read\t$map_read(M)\n";
print "MAP rate\t$map_rate(%)\n";
print "Molecule with PCRdup\t$molecule_with_dup(M)\n";
print "Molecule\t$map_molecule(M)\n";
print "Dup remain\t$dup_rate(%)\n";
close FL;
open FL,"$mpeliup" or die $!;
my $bed_len = 0;
my $bed_1Xcov = 0;
my $bed_100Xcov = 0;
my $bed_500Xcov = 0;
my $bed_1000Xcov = 0;
my $bed_3000Xcov = 0;
my $target_base = 0;
while(<FL>){
	chomp;
	$bed_len++;
	my @arr = split/\t/;
	$target_base += $arr[3];
	next unless($arr[3] >= 1);
	$bed_1Xcov++;
	next unless($arr[3] >= 100);
	$bed_100Xcov++;
	next unless($arr[3] >= 500);
	$bed_500Xcov++;
	next unless($arr[3] >= 1000);
	$bed_1000Xcov++;
	next unless($arr[3] >= 3000);
	$bed_3000Xcov++;
}
close FL;
$bed_1Xcov = sprintf"%.2f",$bed_1Xcov*100/$bed_len;
$bed_100Xcov = sprintf"%.2f",$bed_100Xcov*100/$bed_len;
$bed_500Xcov = sprintf"%.2f",$bed_500Xcov*100/$bed_len;
$bed_1000Xcov = sprintf"%.2f",$bed_1000Xcov*100/$bed_len;
$bed_3000Xcov = sprintf"%.2f",$bed_3000Xcov*100/$bed_len;
my $ave_depth = sprintf"%.2f",$target_base/$bed_len;
my $capture_rate = 0;;
if($outtype eq "CSS"){
        $capture_rate = sprintf"%.2f",$target_base*100/($molecule_with_dup*1000000*($read1_clean_base/$read1_clean_read + $read2_clean_base/$read2_clean_read));
}elsif($outtype eq "ALL"){
        $capture_rate = sprintf"%.2f",$target_base*100/($map_molecule*1000000*($read1_clean_base/$read1_clean_read + $read2_clean_base/$read2_clean_read));
}
$target_base = sprintf"%.2f",$target_base/1000000;
my $data_use = sprintf"%.2f",$target_base*100/$clean_base;
print "Target len\t$bed_len(bp)\n";
print "Target base\t$target_base(M)\n";
print "Capture rate\t$capture_rate(%)\n";
print "Target ave depth\t$ave_depth(X)\n";
print "Target 1X cover\t$bed_1Xcov(%)\n";
print "Target 100X cover\t$bed_100Xcov(%)\n";
print "Target 500X cover\t$bed_500Xcov(%)\n";
print "Target 1000X cover\t$bed_1000Xcov(%)\n";
print "Target 30000X cover\t$bed_3000Xcov(%)\n";
print "Data use\t$data_use(%)\n";

