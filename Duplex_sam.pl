#!/usr/bin/perl -w
use strict;
use Getopt::Long;
=head1 Usage

        perl ConSeqV2.1.pl [option] -i insam -o csssam -O allsam

=head1 Arguments

        -r              Reference genome(default: /thinker/net/ctDNA/hg19/hg19.fa)
        -i              Input sam
        -o              Output CSS sam
        -a              Output ALL sam
        -l              Output stat log
        -e              UMI type(Index/Duplex: default Index)
        -b              UMI length(default 8)
        -c              UMI group allow mismatch(default 1)
        -u              Consensus cutoff(default 0.8 / Set to 0.5 if nonUMI data)
        -s              Sample type(UMI/nonUMI: default UMI)

=cut

my $hg19 = "/thinker/net/ctDNA/hg19/hg19.fa";
my $cutoff = 0.8;
my $sampletype = "UMI";
my $umitype = "Index";
my $umilen = 8;
my $umimis = 1;
my ($insam,$csssam,$allsam,$logfile);

GetOptions(
        "r:s" => \$hg19,
        "i:s" => \$insam,
        "o:s" => \$csssam,
        "a:s" => \$allsam,
        "b:s" => \$umilen,
        "c:s" => \$umimis,
        "e:s" => \$umitype,
        "l:s" => \$logfile,
        "u:s" => \$cutoff,
        "s:s" => \$sampletype,
);

die `pod2text $0` if ($insam eq "" || $csssam eq "" || $allsam eq "");
die `pod2text $0` if ($sampletype ne "UMI" && $sampletype ne "nonUMI");


open FL,"$hg19" or die $!;                           # Read ref genome hash
$/ = ">";
<FL>;
my %chr;
while(<FL>){
        chomp;
        my @arr = split/\n/;
        my $id = shift @arr;
        my $seq = join"",@arr;
        $chr{$id} = $seq;
}
close FL;
$/ = "\n";

my $map_read = 0;
my $dedup_read = 0;
my $dedup_molecule = 0;
my $molecule = 0;
my $molecule_SE_PCR = 0;
my $molecule_SE = 0;
my $molecule_PE_PCR = 0;
my $molecule_PE = 0;
my $mapping_cluster = 0;

open FL,"$insam" or die $!;                  # Read sam file   
open OUT1,">$csssam" or die $!;             # Output css sam file   
open OUT2,">$allsam" or die $!;             # Output all sam file   
my %PEproperly;
my %SEmap;
my %PEunproperly;
my $chrom = "chr1";
while(<FL>){
	chomp;
	if($_ =~ /^@/){
		print OUT1 "$_\n";
		print OUT2 "$_\n";
	}else{
		my @arr = split/\t/;
		my $flag = sprintf"%b",$arr[1];
		next if($flag =~ /1\d{2}$/ || $flag =~ /1\d{8}$/ || $flag =~ /1\d{9}$/ || $arr[5] =~ /H/); # Filter mapping
		$map_read++;
		if($arr[2] ne $chrom){				# End of chromosome
			$mapping_cluster += keys %SEmap;
			$mapping_cluster += keys %PEproperly;
			for my $clu(keys %SEmap){
				&SEprocess($SEmap{$clu});
			}
			for my $clu(keys %PEproperly){
				&PEprocess($PEproperly{$clu});
			}
			%SEmap = ();
			%PEproperly = ();
		}
		my ($end,$cluster);
		if($flag =~ /1\d{3}$/){                                         # SEmap cluster
        	        $end = $arr[3] + &maplen($arr[5]) - 1;
	                $cluster = "$arr[2]:$arr[3]-$arr[2]:$end";
        	        $SEmap{$cluster} .= "$_\n";
		}elsif($flag =~ /1\d{1}$/){					# PE map cluster
			if($flag =~ /0\d{4}$/){
				$end = $arr[3] + abs($arr[8]) - 1;
				$cluster = "$arr[2]:$arr[3]-$arr[2]:$end";
			}else{
				$end = $arr[7] + abs($arr[8]) - 1;
				$cluster = "$arr[2]:$arr[7]-$arr[2]:$end";
			}
			$PEproperly{$cluster} .= "$_\n";
		}elsif($arr[6] =~ /chr/){					# Diff chr cluster
                        if($arr[3] <= $arr[7]){
                                $cluster = "$arr[2]:$arr[3]-$arr[6]:$arr[7]";
                        }else{
                                $cluster = "$arr[6]:$arr[7]-$arr[2]:$arr[3]";
                        }
                        $PEunproperly{$cluster} .= "$_\n";
		}else{								# PE unusual cluster
			
			if($arr[3] <= $arr[7]){
                                $end = $arr[3] + abs($arr[8]) - 1;
                                $cluster = "$arr[2]:$arr[3]-$arr[2]:$end";
                        }else{
                                $end = $arr[7] + abs($arr[8]) - 1;
                                $cluster = "$arr[2]:$arr[7]-$arr[2]:$end";
                        }
			$PEproperly{$cluster} .= "$_\n";
		}
		$chrom = $arr[2];
	}
}
close FL;
$mapping_cluster += keys %SEmap;
$mapping_cluster += keys %PEproperly;
$mapping_cluster += keys %PEunproperly;

for my $clu(keys %SEmap){
	&SEprocess($SEmap{$clu});
}
%SEmap = ();
for my $clu(keys %PEproperly){
        &PEprocess($PEproperly{$clu});
}
%PEproperly = ();
for my $cluster(keys %PEunproperly){
	&PEprocess($PEunproperly{$cluster});
}
%PEunproperly = ();

close OUT1;
close OUT2;

$dedup_read = $molecule_SE * 1 + $molecule_PE * 2;
$dedup_molecule = $molecule_SE + $molecule_PE;
$molecule = $molecule_SE + $molecule_PE;

open OUTL,">$logfile" or die $!;
print OUTL "Mapping read = $map_read\n";
print OUTL "Mapping cluster = $mapping_cluster\n";
print OUTL "Mapping molecule = $molecule\n";
print OUTL "Mapping SE molecule = $molecule_SE\n";
print OUTL "Mapping PE molecule = $molecule_PE\n";
print OUTL "Mapping SE molecule with PCRdup = $molecule_SE_PCR\n";
print OUTL "Mapping PE molecule with PCRdup = $molecule_PE_PCR\n";
print OUTL "Read after dedup = $dedup_read\n";
print OUTL "Molecule after dedup = $dedup_molecule\n";
close OUTL;

sub SEprocess{
	my $cluster = shift @_;
	if($sampletype eq "UMI"){
		my @gop = split/\n/,$cluster;
		my %group;
		if($umitype eq "Duplex"){
			%group = &duplexgroup(@gop);
		}elsif($umitype eq "Index"){
			%group = &indexgroup(@gop);
		}
		for my $k(keys %group){
			$molecule_SE++;
			my $SEread = $group{$k} =~ tr/\n/\n/;
			if($SEread <= 1){
				chomp $group{$k};
				print OUT2 "$group{$k}\tCN:i:1\n";
			}elsif($SEread >=2){
				$molecule_SE_PCR++;
				&dedup($group{$k});
			}
		}
	}elsif($sampletype eq "nonUMI"){
		$molecule_SE++;
		my $SEread = $cluster =~ tr/\n/\n/;
		if($SEread <= 1){
			chomp $cluster;
			print OUT2 "$cluster\tCN:i:1\n";
		}elsif($SEread >=2){
			$molecule_SE_PCR++;
			&dedup($cluster);
		}
	}
}
sub PEprocess{
        my $cluster = shift @_;
        if($sampletype eq "UMI"){
                my @gop = split/\n/,$cluster;
		my %group;
                if($umitype eq "Duplex"){
                        %group = &duplexgroup(@gop);
                }elsif($umitype eq "Index"){
                        %group = &indexgroup(@gop);
                }
                for my $k(keys %group){
			$molecule_PE++;
                        my $PEread = $group{$k} =~ tr/\n/\n/;
                        if($PEread == 2){
				my @crr = split/\n/,$group{$k};
                                print OUT2 "$crr[0]\tCN:i:1\n";
                                print OUT2 "$crr[1]\tCN:i:1\n";
                        }elsif($PEread >=4){
				$molecule_PE_PCR++;
                                &dedup($group{$k});
                        }
                }
        }elsif($sampletype eq "nonUMI"){
		$molecule_PE++;
                my $PEread = $cluster =~ tr/\n/\n/;
                if($PEread == 2){
			my @crr = split/\n/,$cluster;
                        print OUT2 "$crr[0]\tCN:i:1\n";
                        print OUT2 "$crr[1]\tCN:i:1\n";
                }elsif($PEread >= 4){
			$molecule_PE_PCR++;
                        &dedup($cluster);
                }
        }
}

sub duplexgroup{
        my @line = @_;
        my %raw_group;
        my %barcode_read;
        for my $line(@line){
                $line =~ /:UMI_(\w+)_(\w+)/;
		my $barco = "$1"."$2";
                $raw_group{$barco} .= "$line\n";
                $barcode_read{$barco}++;
        }
        for my $aa(sort {$barcode_read{$b} <=> $barcode_read{$a}}keys %barcode_read){
		next unless(exists $barcode_read{$aa});
                delete $barcode_read{$aa};
                for my $bb(keys %barcode_read){
			$raw_group{$bb} =~ /:UMI_(\w+)_(\w+)/;
			my $cc = "$2"."$1";
                        my @xrr = split//,$aa;
                        my @yrr = split//,$bb;
			my @zrr = split//,$cc;
                        my $str_mis1 = 0;
                        for(my $i = 0;$i < $umilen; $i++){
                                if($xrr[$i] ne $yrr[$i]){
                                        $str_mis1++;
                                }
                        }
			my $str_mis2 = 0;
			for(my $i = 0;$i < $umilen; $i++){
                                if($xrr[$i] ne $zrr[$i]){
                                        $str_mis2++;
                                }
                        }
                        if($str_mis1 <= $umimis || $str_mis2 <= $umimis){
                                $raw_group{$aa} .= "$raw_group{$bb}";
                                delete $raw_group{$bb};
                                delete $barcode_read{$bb};
                        }
                }
        }
        return (%raw_group);
}

sub indexgroup{
        my @line = @_;
        my %raw_group;
        my %barcode_read;
        for my $line(@line){
                $line =~ /UMI_(\w+)/;
                $raw_group{$1} .= "$line\n";
                $barcode_read{$1}++;
        }
        for my $aa(sort {$barcode_read{$b} <=> $barcode_read{$a}}keys %barcode_read){
		next unless(exists $barcode_read{$aa});
                delete $barcode_read{$aa};
                for my $bb(keys %barcode_read){
			my @xrr = split//,$aa;
			my @yrr = split//,$bb;
                        my $str_mis = 0;
			for(my $i = 0;$i < $umilen; $i++){
				if($xrr[$i] ne $yrr[$i]){
					$str_mis++;
				}
			}
                        if($str_mis <= $umimis){
                                $raw_group{$aa} .= "$raw_group{$bb}";
                                delete $raw_group{$bb};
                                delete $barcode_read{$bb};
                        }
                }
        }
        return (%raw_group);
}

sub dedup{
	my $cluster = shift @_;
	my @line = split/\n/,$cluster;
	my $read_num = int(($#line + 1)/2);
	my $pre_dedup = "CN:i:$read_num";
	my %postion;
	my %clu_variant;
	my $read_id;
	my %mutpos;
	for my $line(@line){
		my @brr = split/\t/,$line;
		$read_id = $brr[0];
		my $map_end = $brr[3] + &maplen($brr[5]) - 1;
		my @maparr = ($brr[3]..$map_end);
		my %variant = &mutation($brr[5],$brr[12],$brr[3],$brr[9]);
		for my $pp(@maparr){
			$postion{$pp}++;
			if(exists $variant{$pp}){
				$mutpos{$pp}++;
				$clu_variant{$pp}{$variant{$pp}}++;
			}
		}
	}
	my %correct_var;
	for my $pp(keys %mutpos){
			next unless($mutpos{$pp}/$postion{$pp} > $cutoff && $mutpos{$pp} >= 2);
			for my $var(keys %{$clu_variant{$pp}}){
				next unless(exists $clu_variant{$pp}{$var} && $clu_variant{$pp}{$var}/$postion{$pp} > $cutoff && $clu_variant{$pp}{$var} >= 2);
				$correct_var{$pp} = $var;
			}
	}
	for my $line(@line){
		my @brr = split/\t/,$line;
		next unless($read_id eq $brr[0]);
		my $flag = sprintf"%b",$brr[1];
		if($flag =~ /1\d{3}$/){
			$read_num = $#line + 1;
			$pre_dedup = "CN:i:$read_num";
		}
		my ($cigar,$read,$qual,$adj_cigar,$edit_distance,$map_scores) = &rebuild($flag,$brr[2],$brr[3],$brr[5],$brr[9],$brr[10],%correct_var);
		$brr[5] = $cigar;
		$brr[9] = $read;
		$brr[10] = $qual;
		$brr[11] = $edit_distance;
		$brr[12] = $adj_cigar;
		$brr[13] = $map_scores;
		$line = join"\t",@brr;
		print OUT1 "$line\t$pre_dedup\n";
		print OUT2 "$line\t$pre_dedup\n";
	}
}

sub rebuild{
	my ($flag,$chr,$mapstart,$cigar,$read,$qual,%variant) = @_;
	my $map_end = $mapstart + &maplen($cigar) - 1;
	my @maparr = ($mapstart..$map_end);
	my $seq;
	my $seq_qual;
	my $reb_cigar;
	my $reb_adjcigar;
	my $edit_distance = 0;
	my $map_scores = 0;
	
	if($cigar =~ /^(\d+)S/){
		$seq .= substr($read,0,$1);
		$seq_qual .= substr($qual,0,$1);
		$reb_cigar .= "$1S";
		$map_scores -= 5;
	}
	my $step = 0;
	my $step_adj = 0;
	my %skip_pos;
	for my $pos(@maparr){
		next if(exists $skip_pos{$pos});
		if(exists $variant{$pos}){
                        if($variant{$pos} =~ /^(-)(\d+)([A-Z]+)/){
				my @zrr = (1..$2);shift @zrr;           
                                for my $skip(@zrr){                     
                                        my $skip_pos = $pos + $skip - 1;
                                        $skip_pos{$skip_pos}++;         
                                }
                                $edit_distance += $2;
                                $map_scores -= (6 + $2);
                                $reb_cigar .= "$step"."M"."$2"."D";
                                $step = 0;
                                $reb_adjcigar .= "$step_adj"."^"."$3";
                                $step_adj = 0;
                        }elsif($variant{$pos} =~ /^(\+)(\d+)([A-Z]+)/){
                                $edit_distance += $2;
                                $map_scores -= (6 + $2);
                                my $chr_pos = substr($chr{$chr},$pos - 1,1);
                                $seq .= "$chr_pos"."$3";
                                my $ins_base = $3;
                                $ins_base =~ tr/[ATCGN]/[IIIII]/;
                                $seq_qual .= "I"."$ins_base";
                                $step++;
                                $reb_cigar .= "$step"."M"."$2"."I";
                                $step = 0;
                                $step_adj++;
                        }elsif($variant{$pos} =~ /([A-Z]+)>([A-Z]+)/){
                                $edit_distance++;
                                $map_scores -= 4;
                                $seq .= "$2";
                                $seq_qual .= "I";
                                $step++;
                                $reb_adjcigar .= "$step_adj"."$1";
                                $step_adj = 0;
                        }
                }else{
			$map_scores++;
                        $step++;
                        $step_adj++;
                        $seq .= substr($chr{$chr},$pos - 1,1);
                        $seq_qual .= "I";
                }
	}
	$reb_cigar .= "$step"."M";
	$reb_adjcigar .= "$step_adj";
	$reb_adjcigar = "MD:Z:$reb_adjcigar";
	if($cigar =~ /(\d+)S$/){
		my $sta = length($read) - $1;
		$seq .= substr($read,$sta,$1);
		$seq_qual .= substr($qual,$sta,$1);
		$reb_cigar .= "$1S";
		$map_scores -= 5;
	}
	$seq = uc($seq);
	$edit_distance = "NM:i:$edit_distance";
	if($flag =~ /1\d{3}$/){
		$map_scores -= 17;
	}
	$map_scores = "AS:i:$map_scores";
	return ($reb_cigar,$seq,$seq_qual,$reb_adjcigar,$edit_distance,$map_scores);
}

sub mutation{
	my ($cigar,$adj_cigar,$mapstart,$read) = @_;
	my %variant;
	my $dist = $mapstart - 1;
	while($cigar =~ /(\d+)(I|D|M|S)/g){
                my $extend = $1;
                my $cigar_half = $`;
                my $mut_type = $2;
                if($mut_type eq "I"){
			my $len_tmp = 0;
			while($cigar_half =~ /(\d+)(M|S|I)/g){
				$len_tmp += $1;
			}
			my $ins_base = substr($read,$len_tmp,$extend);
                        $variant{$dist} = "+$extend"."$ins_base";
                }elsif($mut_type eq "D"){
                        my $del_dist = $dist + 1;
                        $variant{$del_dist} = "-$extend";
                        $dist += $extend;
                }elsif($mut_type eq "M"){
                        $dist += $extend;
                }
        }
        my $dist_adj = $mapstart - 1;
        while($adj_cigar =~ /(\d+)(A|T|C|G|\^)/g){
		my $cigar_half = $';
		my $extend = $1;
		my $mut_type = $2;
                if($mut_type eq "^"){
                        $dist_adj += ($extend + 1);
			$cigar_half =~ /^([ATCG]+)/;
			my $del_seq = $1;
			$variant{$dist_adj} .= $del_seq;
			my $del_len = length $del_seq;
                        $dist_adj += ($del_len - 1);
                }else{
                        $dist_adj += ($extend + 1);
			my $mis_dist = $dist_adj - $mapstart;
			my $mis_base = &mismatch($mis_dist,$cigar,$read);
                        $variant{$dist_adj} = "$mut_type".">"."$mis_base";
                }
        }
	return %variant;
}

sub mismatch{
	my ($mic_dist,$cigar,$read) = @_;
	my $del_len = 0;
	my $ins_len = 0;
	my $clip_len = 0;
	my $dis_read = 0;
	while($cigar =~ /(\d+)(S|I|D|M)/g){
		next if($dis_read > $mic_dist);
		if($2 eq "S"){
			$clip_len += $1;
		}elsif($2 eq "I"){
			$ins_len += $1;
		}elsif($2 eq "D"){
			$dis_read += $1;
			$del_len += $1;
		}elsif($2 eq "M"){
			$dis_read += $1;
		}
	}
	$mic_dist = $mic_dist + $clip_len + $ins_len - $del_len;
	my $mis_base = substr($read,$mic_dist,1);
	return $mis_base;
}

sub maplen{
        my $cigar = shift @_;
        my $len = 0;
        while($cigar =~ /(\d+)(M|D)/g){
                $len += $1;
        }
        return ($len);
}

