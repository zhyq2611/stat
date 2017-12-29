#!/usr/bin/perl -w
use strict;
use Getopt::Long;
=head1 Usage

        perl Duplex_sam.pl [option] <bwamem.sam> <out.sam>

=head1 Arguments

        -r              Reference genome(default /thinker/net/ctDNA/hg19/hg19.fa)
        -i              Bwa mem out put sam
        -o              Out put sam
        -t              OUT put read type: only Consensus or all read(CSS or ALL: default CSS)
        -u              Consensus cutoff: threshold for a variant in pcr dup clster or set to WT(default 0.5)
        -s              Sample type UMI or nonUMI(for nonUMI sample build Consensus sequence base on map start and end)

=cut

my $hg19 = "/thinker/net/ctDNA/hg19/hg19.fa";
my $type = "CSS";
my $cutoff = 0.5;
my $sampletype = "UMI";
my ($sam,$outsam);

GetOptions(
        "r:s" => \$hg19,
        "i:s" => \$sam,
        "o:s" => \$outsam,
        "t:s" => \$type,
        "u:s" => \$cutoff,
        "s:s" => \$sampletype,
);
die `pod2text $0` if ($sam eq "" || $outsam eq "");

open FL,"$hg19" or die $!;
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

my $map_Pe_pro = 0;
my $map_Se = 0;
my $map_diff_chr = 0;
my $map_Pe_unpro = 0;

open FL,"$sam" or die $!;
open OUT,">$outsam" or die $!;
my %PEproperly;
my %SEmap;
my %PEunproperly;
while(<FL>){
	chomp;
	if($_ =~ /^@/){
		print OUT "$_\n";
	}else{
		my @arr = split/\t/;
		my $flag = sprintf"%b",$arr[1];
		next if($flag =~ /0$/ || $flag =~ /1\d{2}$/ || $flag =~ /1\d{8}$/ || $flag =~ /1\d{9}$/ || $arr[5] =~ /H/);
		my $cluster;
		if($flag =~ /1\d{3}$/){
			$map_Se++;
        	        my $end = $arr[3] + &maplen($arr[5]) - 1;
	                $cluster = "$arr[2]:$arr[3]-$arr[2]:$end";
        	        $SEmap{$cluster} .= "$_\n";
	        }elsif($arr[6] =~ /chr/ || $flag =~ /00\d{4}$/ || $flag =~ /11\d{4}$/){
			if($arr[6] eq "="){
                	        $arr[6] = $arr[2];
				$map_Pe_unpro++;
	                }else{
				$map_diff_chr++;
			}       
        	        if($arr[3] <= $arr[7]){
                	        $cluster = "$arr[2]:$arr[3]-$arr[6]:$arr[7]";
	                }else{  
        	                $cluster = "$arr[6]:$arr[7]-$arr[2]:$arr[3]";
                	}
	                $PEunproperly{$cluster} .= "$_\n";
		}else{
			$map_Pe_pro++;
			if($flag =~ /0\d{4}$/){
				my $end = $arr[3] + abs($arr[8]) - 1;
				$cluster = "$arr[2]:$arr[3]-$arr[2]:$end";
				$PEproperly{$cluster} .= "$_\n";
			}else{
				my $end = $arr[7] + abs($arr[8]) - 1;
				$cluster = "$arr[2]:$arr[7]-$arr[2]:$end";
				$PEproperly{$cluster} .= "$_\n";
			}	
		}
	}
}
close FL;

if($sampletype eq "UMI"){
	my $molecule_nopcr = 0;
	my $molecule_pcr = 0;
	my $map_clu = 0;
	my $map_clu_single = 0;
	my $map_clu_multi = 0;
	my $map_clu_multi_barcode = 0;
	for my $cluster(keys %PEproperly){
		$map_clu++;
		my $line = $PEproperly{$cluster} =~ tr/\n/\n/;
		if($line == 2){
			$map_clu_single++;
			$molecule_nopcr++;
			next unless($type eq "ALL");
			my @crr = split/\n/,$PEproperly{$cluster};
			print OUT "$crr[0]\tCN:i:1\n";
			print OUT "$crr[1]\tCN:i:1\n";
		}else{
			my @gop = split/\n/,$PEproperly{$cluster};
			my %group =  &group(@gop);
			my $bar_clu = keys %group;
			if($bar_clu == 1){
				$map_clu_single++;
			}else{
				$map_clu_multi++;
				$map_clu_multi_barcode += $bar_clu;
			}
			for my $k(keys %group){
				my $PEread = $group{$k} =~ tr/\n/\n/;
				if($PEread == 2){
					$molecule_nopcr++;
					next unless($type eq "ALL");
					my @crr = split/\n/,$group{$k};
					print OUT "$crr[0]\tCN:i:1\n";
	                                print OUT "$crr[1]\tCN:i:1\n";
				}else{
					$molecule_pcr++;
                        	        &dedup($group{$k});
	                        }
			}
		}
	
	}
	for my $cluster(keys %PEunproperly){
		$map_clu++;
		my $line = $PEunproperly{$cluster} =~ tr/\n/\n/;
		if($line == 2){
			$map_clu_single++;
			$molecule_nopcr++;
			next unless($type eq "ALL");
			my @crr = split/\n/,$PEunproperly{$cluster};
			print OUT "$crr[0]\tCN:i:1\n";
			print OUT "$crr[1]\tCN:i:1\n";
	        }else{
			my @gop = split/\n/,$PEunproperly{$cluster};
			my %group =  &group(@gop);
			my $bar_clu = keys %group;
	                if($bar_clu == 1){
	                        $map_clu_single++;
	                }else{
	                        $map_clu_multi++;
	                        $map_clu_multi_barcode += $bar_clu;
	                }
			for my $k(keys %group){
				my $PEread = $group{$k} =~ tr/\n/\n/;
				if($PEread == 2){
					$molecule_nopcr++;
					next unless($type eq "ALL");
					my @crr = split/\n/,$group{$k};
					print OUT "$crr[0]\tCN:i:1\n";
					print OUT "$crr[1]\tCN:i:1\n";
				}else{
					$molecule_pcr++;
					&dedup($group{$k});
				}
			}
	        }
	}
	for my $cluster(keys %SEmap){
		$map_clu++;
		my $line = $SEmap{$cluster} =~ tr/\n/\n/;
		if($line == 1){
			$map_clu_single++;
			$molecule_nopcr++;
			next unless($type eq "ALL");
			$SEmap{$cluster} =~ s/\n/$2/;
	                print OUT "$SEmap{$cluster}\tCN:i:1\n";
	        }else{
			my @gop = split/\n/,$SEmap{$cluster};
	                my %group =  &group(@gop);
			my $bar_clu = keys %group;
	                if($bar_clu == 1){
	                        $map_clu_single++;
	                }else{
	                        $map_clu_multi++;
	                        $map_clu_multi_barcode += $bar_clu;
	                }
	                for my $k(keys %group){
	                        my $PEread = $group{$k} =~ tr/\n/\n/;
	                        if($PEread == 1){
					$molecule_nopcr++;
					next unless($type eq "ALL");
					$group{$k} =~ s/\n//;
	                                print OUT "$group{$k}\tCN:i:1\n";
	                        }else{
					$molecule_pcr++;
	                                &dedup($group{$k});
	                        }
	                }
	        }
	}
	print "Read map in Pair and Diff Chr\t$map_diff_chr\n";
	print "Read map in pair properly\t$map_Pe_pro\n";
	print "Read map in pair unproperly\t$map_Pe_unpro\n";
	print "Read map in single\t$map_Se\n";
	print "Molecule with dup\t$molecule_pcr\n";
	print "Molecule without dup\t$molecule_nopcr\n";
	print "Map cluster\t$map_clu\n";
	print "Map cluster single Molecule\t$map_clu_single\n";
	print "Map cluster multi Molecule\t$map_clu_multi\n";
	print "Molecule in Map cluster multi Molecule\t$map_clu_multi_barcode\n";
}elsif($sampletype eq "nonUMI"){
	my $molecule_nopcr = 0;
	my $molecule_pcr = 0;
	for my $cluster(keys %PEproperly){
	        my $line = $PEproperly{$cluster} =~ tr/\n/\n/;
	        if($line == 2){
	                $molecule_nopcr++;
			next unless($type eq "ALL");
	                my @crr = split/\n/,$PEproperly{$cluster};
	                print OUT "$crr[0]\tCN:i:1\n";
	                print OUT "$crr[1]\tCN:i:1\n";
	        }else{
        	        $molecule_pcr++;
	                &dedup($PEproperly{$cluster});
	        }
	}
	for my $cluster(keys %PEunproperly){
	        my $line = $PEunproperly{$cluster} =~ tr/\n/\n/;
	        if($line == 2){
	                $molecule_nopcr++;
			next unless($type eq "ALL");
	                my @crr = split/\n/,$PEunproperly{$cluster};
	                print OUT "$crr[0]\tCN:i:1\n";
        	        print OUT "$crr[1]\tCN:i:1\n";
	        }else{
	                $molecule_pcr++;
	                &dedup($PEunproperly{$cluster});
	        }
	}
	for my $cluster(keys %SEmap){
	        my $line = $SEmap{$cluster} =~ tr/\n/\n/;
	        if($line == 1){
	                $molecule_nopcr++;
			next unless($type eq "ALL");
	                print OUT "$SEmap{$cluster}\tCN:i:1\n";
	        }else{
	                $molecule_pcr++;
	                &dedup($SEmap{$cluster});
	        }
	}
	print "Read map in Pair and Diff Chr\t$map_diff_chr\n";
        print "Read map in pair properly\t$map_Pe_pro\n";
        print "Read map in pair unproperly\t$map_Pe_unpro\n";
        print "Read map in single\t$map_Se\n";
        print "Molecule with dup\t$molecule_pcr\n";
        print "Molecule without dup\t$molecule_nopcr\n";
}
close OUT;


sub group{
        my @line = @_;
        my %raw_group;
        my %barcode_read;
        for my $line(@line){
                $line =~ /UMI_(\w+)/;
		my $barcode = $1;
                $raw_group{$barcode} .= "$line\n";
                $barcode_read{$barcode}++;
        }
        for my $aa(sort {$barcode_read{$b} <=> $barcode_read{$a}}keys %barcode_read){
                delete $barcode_read{$aa};
                for my $bb(sort {$barcode_read{$b} <=> $barcode_read{$a}}keys %barcode_read){
                        my $str_match = &str_match($aa,$bb);
                        if($str_match >= 7){
                                $raw_group{$aa} .= "$raw_group{$bb}";
                                delete $raw_group{$bb};
                                delete $barcode_read{$bb};
                        }
                }
        }
        return (%raw_group);
}
sub str_match{
        my ($str1,$str2) = @_;
        my @xrr = split//,$str1;
        my @yrr = split//,$str2;
        my $str = 0;
        my $str_match = 0;
        for my $nn(@xrr){
                if($nn eq $yrr[$str]){
                        $str_match++;
                }
                $str++;
        }
        return($str_match);
}

sub dedup{
	my $clutser = shift @_;
	chomp $clutser;
	my @line = split/\n/,$clutser;
	my $read_num = int(($#line + 1)/2);
	my $pre_dedup = "CN:i:$read_num";
	my %postion;
	my %clu_variant;
	my $read_id;
	for my $line(@line){
		my @brr = split/\t/,$line;
		$read_id = $brr[0];
		my $map_end = $brr[3] + &maplen($brr[5]) - 1;
		my @maparr = ($brr[3]..$map_end);
		my %variant = &mutation($brr[5],$brr[12],$brr[3],$brr[9]);
		for my $pp(@maparr){
			$postion{$pp}++;
			if(exists $variant{$pp}){
				$clu_variant{$pp}{$variant{$pp}}++;
			}
		}
	}
	my %correct_var;
	for my $pp(keys %postion){
			for my $var(keys %{$clu_variant{$pp}}){
				next unless(exists $clu_variant{$pp}{$var} && $clu_variant{$pp}{$var}/$postion{$pp} > $cutoff && $clu_variant{$pp}{$var} >= 2);
				$correct_var{$pp} = $var;
			}
	}
	for my $line(@line){
		my @brr = split/\t/,$line;
		my $flag = sprintf"%b",$brr[1];
		next unless($read_id eq $brr[0]);
		my $map_end = $brr[3] + &maplen($brr[5]) - 1;
		my @maparr = ($brr[3]..$map_end);
		my ($cigar,$read,$qual,$adj_cigar,$edit_distance,$map_scores) = &rebuild($flag,$brr[2],$brr[3],$brr[5],$brr[9],$brr[10],%correct_var);
		$brr[5] = $cigar;
		$brr[9] = $read;
		$brr[10] = $qual;
		$brr[11] = $edit_distance;
		$brr[12] = $adj_cigar;
		$brr[13] = $map_scores;
		$line = join"\t",@brr;
		print OUT "$line\t$pre_dedup\n";
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
                my $mut_type = $2;
                if($mut_type eq "I"){
			my $ins_base = substr($read,&inspos($`),$extend);
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
                if($2 eq "^"){
                        $dist_adj += ($1 + 1);
			$variant{$dist_adj} .= &delpos($');
			my $del_len = length &delpos($');
                        $dist_adj += ($del_len - 1);
                }else{
                        $dist_adj += ($1 + 1);
			my $mis_dist = $dist_adj - $mapstart;
			my $mis_base = &mismatch($mis_dist,$cigar,$read);
                        $variant{$dist_adj} = "$2".">"."$mis_base";
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

sub inspos{
	my $cigar = shift @_;
	my $len = 0;
	while($cigar =~ /(\d+)(M|S|I)/g){
                $len += $1;
        }
        return ($len);
}	
sub delpos{
	my $cigar = shift @_;
	$cigar =~ /^([ATCG]+)/;
	return $1;
}

sub maplen{
        my $cigar = shift @_;
        my $len = 0;
        while($cigar =~ /(\d+)(M|D)/g){
                $len += $1;
        }
        return ($len);
}


