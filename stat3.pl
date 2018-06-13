#!/usr/bin/perl -w
use strict;
use Getopt::Long;
=head1 Usage

        perl stat.pl -s snp -i indel -r refscan -a altscan -f fuisonlog

=head1 Arguments

        -s              Snp file(WS-P1_css.snp.txt)
        -i              Indel file(WS-P1_css.indel.txt)
        -r		Mutscan ref file(WS-P1_mutscan_ref.html)
        -a              Mutscan alt file(WS-P1_mutscan_alt.html)
        -f              Genefuse fusion log file(WS-P1.fusion.log)

=cut
my ($snp,$indel,$refscan,$altscan,$fusion);

GetOptions(
        "s:s" => \$snp,
        "i:s" => \$indel,
        "r:s" => \$refscan,
        "a:s" => \$altscan,
        "f:s" => \$fusion,
);

die `pod2text $0` if ($snp eq "" || $indel eq "" || $refscan eq "" || $altscan eq "" || $fusion eq "");

my $L858R = "Negative";
my $T790M = "Negative";
my $del19 = "Negative";
my $ALK = "Negative";
my $L858R_ref = 0;
my $T790M_ref = 0;
my $del19_ref = 0;
my $ALK_ref_m = 0;
my $ALK_ref_f = 0;
my $L858R_alt = 0;
my $T790M_alt = 0;
my $del19_alt = 0;
my $ALK_alt_m = 0;
my $ALK_alt_f = 0;

my $L858R_depth = 0;
my $T790M_depth = 0;
my $del19_depth = 0;
my $ALK_depth = 0;


my $L858R_depth_ref = 0;
my $T790M_depth_ref = 0;
my $del19_depth_ref = 0;
my $ALK_depth_ref = 0;
my $ALK_depth_ref_m = 0;
my $ALK_depth_ref_f = 0;

my $L858R_depth_alt = 0;
my $T790M_depth_alt = 0;
my $del19_depth_alt = 0;
my $ALK_depth_alt = 0;
my $ALK_depth_alt_m = 0;
my $ALK_depth_alt_f = 0;


open FL,"$snp" or die $!;
my $title = <FL>;
my @ti = split/\t/,$title;
while(<FL>){
	chomp;
	my %hash;
	@hash{@ti} = my @arr = split/\t/;
	my $overlap = $hash{Multiple_Overlap_alt} + $hash{One_Overlap_alt};
	my $molecule = $hash{Multiple_Overlap_alt} + $hash{One_Overlap_alt} + $hash{Multiple_Single_alt} + $hash{One_Single_alt};
	next unless($overlap >= 2 || ($molecule >= 3 && $overlap >= 1));
	$hash{Tumor_AF} =~ s/%//;
	next unless($hash{Tumor_AF} >= 0.2);
	if($hash{Gene} eq "EGFR" && $hash{Protein} =~ /p\.T790M/){
			$T790M = "Positive";
			$T790M_alt = $molecule;
			$T790M_depth_alt = $hash{Tumor_altread};
	}elsif($hash{Gene} eq "EGFR" && $hash{Protein} =~ /p\.L858R/){
			$L858R = "Positive";
			$L858R_alt = $molecule;
			$L858R_depth_alt = $hash{Tumor_altread};
	}
}
close FL;

open FL,"$indel" or die $!;
<FL>;
while(<FL>){
        chomp;
        my %hash;
        @hash{@ti} = my @arr = split/\t/;
        my $overlap = $hash{Multiple_Overlap_alt} + $hash{One_Overlap_alt};
        my $molecule = $hash{Multiple_Overlap_alt} + $hash{One_Overlap_alt} + $hash{Multiple_Single_alt} + $hash{One_Single_alt};
	next unless($overlap >= 2 || ($molecule >= 3 && $overlap >= 1));
	$hash{Tumor_AF} =~ s/%//;
	next unless($hash{Tumor_AF} >= 0.2);
        if($hash{Gene} eq "EGFR" && $hash{Nucleotide} =~ /c\.2235_2249del/){
		$del19 = "Positive";
		$del19_alt = $molecule;
		$del19_depth_alt = $hash{Tumor_altread};
        }
}
close FL;
open FL,"$fusion" or die $!;
while(<FL>){
	chomp;
	next unless($_ =~ /EML4_ENST00000318522\.5:intron:13(.*)ALK_ENST00000389048\.3:intron:19(.*)\(total: (\d+), unique:(\d+)\)/);
	$ALK_alt_f = $4;
	$ALK_depth_alt_f = $3;	
}
close FL;


open FL,"$refscan" or die $!;
while(<FL>){
	chomp;
	if($_ =~ /EGFR_L858R \((\d+) reads support, (\d+) unique\)/){
		$L858R_ref = $2;
		$L858R_depth_ref = $1;
        }elsif($_ =~ /EGFR_T790M \((\d+) reads support, (\d+) unique\)/){
                $T790M_ref = $2;
		$T790M_depth_ref = $1;
        }elsif($_ =~ /EGFR_19del \((\d+) reads support, (\d+) unique\)/){
                $del19_ref = $2;
		$del19_depth_ref = $1;
        }elsif($_ =~ /EML4_29447649_ALK \((\d+) reads support, (\d+) unique\)/){
                $ALK_ref_m = $2;
		$ALK_depth_ref_m = $1;
        }elsif($_ =~ /EML4_29447818_ALK \((\d+) reads support, (\d+) unique\)/){
		$ALK_ref_f = $2;
		$ALK_depth_ref_f = $1;
	}

}
close FL;

my $EML4_ALK_start = 0;
my $EML4_ALK_start_depth = 0;
my $EML4_ALK_end = 0;
my $EML4_ALK_end_depth = 0;
open FL,"$altscan" or die $!;
while(<FL>){
	chomp;
        if($_ =~ /EGFR_L858R \((\d+) reads support, (\d+) unique\)/){
		if($L858R eq "Positive" && $L858R_alt < $2){
			$L858R_alt = $2;
			$L858R_depth_alt = $1;
		}
        }elsif($_ =~ /EGFR_T790M \((\d+) reads support, (\d+) unique\)/){
		if($T790M eq "Positive" && $T790M_alt < $2){
                        $T790M_alt = $2;
			$T790M_depth_alt = $1;
                }
        }elsif($_ =~ /EGFR_19del \((\d+) reads support, (\d+) unique\)/){
		if($del19 eq "Positive" && $del19_alt < $2){
                        $del19_alt = $2;
			$del19_depth_alt = $1;
                }
        }elsif($_ =~ /EML4_ALK_start \((\d+) reads support, (\d+) unique\)/){
		$EML4_ALK_start = $2;	
		$EML4_ALK_start_depth = $1;
        }elsif($_ =~ /EML4_ALK_end \((\d+) reads support, (\d+) unique\)/){
		$EML4_ALK_end = $2;
		$EML4_ALK_end_depth = $1;
	}
}
close FL;
$ALK_alt_m = $EML4_ALK_start + $EML4_ALK_end;
$ALK_depth_alt_m = $EML4_ALK_end_depth + $EML4_ALK_start_depth;

my $ALK_alt = 0;
my $ALK_ref = 0;
if($ALK_alt_f >= $ALK_alt_m){
	$ALK_alt = $ALK_alt_f;
	$ALK_ref = $ALK_ref_f;
	$ALK_depth_alt = $ALK_depth_alt_f;
	$ALK_depth_ref = $ALK_depth_ref_f;
}else{
	$ALK_alt = $ALK_alt_m;
	$ALK_ref = $ALK_ref_m;
	$ALK_depth_alt = $ALK_depth_alt_m;
	$ALK_depth_ref = $ALK_depth_ref_m;
}
if($ALK_alt >= 3 && $ALK_depth_alt > $ALK_alt){
	$ALK = "Positive";
}else{
	$ALK_alt = 0;
	$ALK_depth_alt = 0;
}

my $L858R_VAF = sprintf"%.2f",$L858R_alt*100/($L858R_alt + $L858R_ref);
my $T790M_VAF = sprintf"%.2f",$T790M_alt*100/($T790M_alt + $T790M_ref);
my $del19_VAF = sprintf"%.2f",$del19_alt*100/($del19_alt + $del19_ref);
my $ALK_VAF = sprintf"%.2f",$ALK_alt*100/($ALK_alt + $ALK_ref);

my $L858R_d = $L858R_alt + $L858R_ref;
my $T790M_d = $T790M_alt + $T790M_ref;
my $del19_d = $del19_alt + $del19_ref;
my $ALK_d = $ALK_alt + $ALK_ref;

$L858R_depth = $L858R_depth_ref + $L858R_depth_alt;
$T790M_depth = $T790M_depth_ref + $T790M_depth_alt;
$del19_depth = $del19_depth_ref + $del19_depth_alt;
$ALK_depth = $ALK_depth_ref + $ALK_depth_alt;

print "基因\t突变\t检测结果\t突变DNA分子\t突变比例\t总DNA分子\t测序深度(X)\n";
print "EGFR\tL858R\t$L858R\t$L858R_alt\t$L858R_VAF\t$L858R_d\t$L858R_depth\n";
print "EGFR\tT790M\t$T790M\t$T790M_alt\t$T790M_VAF\t$T790M_d\t$T790M_depth\n";
print "EGFR\t19del\t$del19\t$del19_alt\t$del19_VAF\t$del19_d\t$del19_depth\n";
print "ALK\tEML4-ALK\t$ALK\t$ALK_alt\t$ALK_VAF\t$ALK_d\t$ALK_depth\n";
