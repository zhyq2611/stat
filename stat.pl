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
my $ALK_ref = 0;
my $L858R_alt = 0;
my $T790M_alt = 0;
my $del19_alt = 0;
my $ALK_alt = 0;
my $L858R_molecule = 0;
my $T790M_molecule = 0;
my $del19_molecule = 0;
my $ALK_molecule = 0;

open FL,"$snp" or die $!;
my $title = <FL>;
my @ti = split/\t/,$title;
while(<FL>){
	chomp;
	my %hash;
	@hash{@ti} = my @arr = split/\t/;
	my $overlap = $hash{Multiple_Overlap_alt} + $hash{One_Overlap_alt};
	my $molecule = $hash{Multiple_Overlap_alt} + $hash{One_Overlap_alt} + $hash{Multiple_Single_alt} + $hash{One_Single_alt};
	next unless($overlap >= 1 && $molecule >= 2);
	if($hash{Gene} eq "EGFR" && $hash{Protein} =~ /p\.T790M/){
			$T790M = "Positive";
	}elsif($hash{Gene} eq "EGFR" && $hash{Protein} =~ /p\.L858R/){
			$L858R = "Positive";
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
	next unless($overlap >= 1 && $molecule >= 2);
        if($hash{Gene} eq "EGFR" && $hash{Nucleotide} =~ /c\.2235_2249del/){
		$del19 = "Positive";
        }
}
close FL;
open FL,"$fusion" or die $!;
while(<FL>){
	chomp;
	next unless($_ =~ /EML4_ENST00000318522\.5:intron:13(.*)ALK_ENST00000389048\.3:intron:19(.*)\(total: (\d+), unique:(\d+)\)/);
	next unless($4 >= 2);
	$ALK = "Positive";
	$ALK_alt = $3;
	$ALK_molecule = $4;
}
close FL;


open FL,"$refscan" or die $!;
while(<FL>){
	chomp;
	if($_ =~ /EGFR_L858R \((\d+) reads support, (\d+) unique\)/){
		$L858R_ref = $1;
        }elsif($_ =~ /EGFR_T790M \((\d+) reads support, (\d+) unique\)/){
                $T790M_ref = $1;
        }elsif($_ =~ /EGFR_19del \((\d+) reads support, (\d+) unique\)/){
                $del19_ref = $1;
        }elsif($_ =~ /EML4_ALK \((\d+) reads support, (\d+) unique\)/){
                $ALK_ref = $1;
        }

}
close FL;

open FL,"$altscan";
while(<FL>){
	chomp;
        if($_ =~ /EGFR_L858R \((\d+) reads support, (\d+) unique\)/){
		if($L858R eq "Positive"){
			$L858R_alt = $1;
			$L858R_molecule = $2;
		}elsif($2 >= 3 && $1 > $2){	
			$L858R = "Positive";
			$L858R_alt = $1;
			$L858R_molecule = $2;
		}
        }elsif($_ =~ /EGFR_T790M \((\d+) reads support, (\d+) unique\)/){
		if($T790M eq "Positive"){
                        $T790M_alt = $1;
			$T790M_molecule = $2;
                }elsif($2 >= 3 && $1 > $2){
			$T790M = "Positive";
                        $T790M_alt = $1;
			$T790M_molecule = $2;
                }
        }elsif($_ =~ /EGFR_19del \((\d+) reads support, (\d+) unique\)/){
		if($del19 eq "Positive"){
                        $del19_alt = $1;
			$del19_molecule = $2;
                }elsif($2 >= 3 && $1 > $2){
			$del19 = "Positive";
                        $del19_alt = $1;
			$del19_molecule = $2;
                }
        }elsif($_ =~ /EML4_ALK_start \((\d+) reads support, (\d+) unique\)/){
		$ALK = "Positive";
		$ALK_alt += $1;	
		$ALK_molecule += $2;
        }elsif($_ =~ /EML4_ALK_end \((\d+) reads support, (\d+) unique\)/){
		$ALK = "Positive";
		$ALK_alt += $1;
		$ALK_molecule += $2;
	}
}
close FL;
my $L858R_VAF = sprintf"%.2f",$L858R_alt*100/($L858R_alt + $L858R_ref);
my $T790M_VAF = sprintf"%.2f",$T790M_alt*100/($T790M_alt + $T790M_ref);
my $del19_VAF = sprintf"%.2f",$del19_alt*100/($del19_alt + $del19_ref);
my $ALK_VAF = sprintf"%.2f",$ALK_alt*100/($ALK_alt + $ALK_ref);

print "基因\t突变\t检测结果\t突变DNA分子\t突变比例\t突变read\t野生型read\n";
print "EGFR\tL858R\t$L858R\t$L858R_molecule\t$L858R_VAF\t$L858R_alt\t$L858R_ref\n";
print "EGFR\tT790M\t$T790M\t$T790M_molecule\t$T790M_VAF\t$T790M_alt\t$T790M_ref\n";
print "EGFR\t19del\t$del19\t$del19_molecule\t$del19_VAF\t$del19_alt\t$del19_ref\n";
print "ALK\tEML4-ALK\t$ALK\t$ALK_molecule\t$ALK_VAF\t$ALK_alt\t$ALK_ref\n";
