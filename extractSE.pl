#!/usr/bin/perl -w
use strict;
my ($infile,$outfile) = @ARGV;
=head1 Usage

        perl extract.pl  <infile>  <outfile>

=head1 Arguments

        <infile>        Single sample MrBam result
        <outfile>       out put file

=cut

die `pod2text $0` if (@ARGV != 2);

open FL,"/thinker/net/ctDNA/Gene_transcript.list" or die $!;
my %genetran;
my %tran_info;
while(<FL>){
        chomp;
        my @arr = split/\t/;
        $genetran{$arr[0]} = $arr[1];
        $tran_info{$arr[0]} = $arr[2];
}
close FL;

open FL,"$infile" or die $!;
open FLS,">$outfile" or die $!;
my $title = <FL>;
print FLS "Chr\tPosition\tRef\tAlt\tdbSNP\tGene\tTranscript\tRegion\tNucleotide\tProtein\tSIFT\tClinVar\tCOSMIC_ID\tCOSMIC_OCCU\tESP6500_AF\tKgenome_AF\tTumor_AF\tTumor_altread\tMultiple_Overlap_ref\tMultiple_Single_ref\tOne_Overlap_ref\tOne_Single_ref\tMultiple_Overlap_alt\tMultiple_Single_alt\tOne_Overlap_alt\tOne_Single_alt\tTumor_depth\tCOUNT\tMEAN_SUPPORT\tMEAN_FREQUENCY\n";
while(<FL>){
	chomp;
	my @arr = split/\t/;
	my $MEAN_FREQUENCY = pop @arr;
	my $MEAN_SUPPORT = pop @arr;
	my $COUNT = pop @arr;
	my $Mut = "$arr[0]\t$arr[1]\t$arr[3]\t$arr[4]";
	my @zrr = split/\:/,$arr[-1];
	my @xrr = split/\,/,$zrr[-1];
	my $tumor_depth = $zrr[3];
	my $tumor_altread = $zrr[5];
	my $tumor_AF = $zrr[6];
	my $Multiple_Overlap_ref = $xrr[0];
        my $Multiple_Single_ref = $xrr[1] + $xrr[2];
        my $One_Overlap_ref = $xrr[3];
        my $One_Single_ref = $xrr[4] + $xrr[5];
        my $Multiple_Overlap_alt = $xrr[6];
        my $Multiple_Single_alt = $xrr[7] + $xrr[8];
        my $One_Overlap_alt = $xrr[9];
        my $One_Single_alt = $xrr[10] + $xrr[11];

	my $EPS6500_AF = $arr[12];
	my $kgenome_AF = $arr[13];
	my $dbSNP = $arr[17];
	my $SIFT = $arr[18];
	my $clinvar = $arr[44];

	my $cosmid = "NA";
	my $cosmic_occu = 0;
        if($arr[43] =~ /ID=COSM/){
                my @crr = split/;/,$arr[43];
                $crr[0] =~ s/ID=//;
                $cosmid = $crr[0];
		$crr[1] =~ s/OCCURENCE=//;
		$crr[1] =~ s/\(.*?\)//g;
		my @drr = split/\,/,$crr[1];
		for my $nn(@drr){
        	        $nn =~ s/\)//;
	                $cosmic_occu += $nn;
	        }
	}
	my $annotation = "NA";
	if($arr[8] =~ /unknown/i || $arr[9] =~ /wholegene/){
                $annotation = "$arr[6]\t.\t.\t.\tnil";
        }else{
		if($arr[5] eq "splicing" && $arr[7] ne "."){
	                my @crr = split/,|;/,$arr[7];
	                my ($tran,$exon,$nt_change) = split/:/,$crr[0];
	                for my $anno_line(@crr){
	                        my @drr = split/\:/,$anno_line;
	                        next unless(exists $genetran{$arr[6]} && $genetran{$arr[6]} eq $drr[0]);
	                        ($tran,$exon,$nt_change) = split/\:/,$anno_line;
	                }
	                $exon =~ s/exon//;
	                my $intron = "unknown";
	                if($nt_change =~ /c\.(\d+)(\-|\+)(\d+)/ && exists $tran_info{$arr[6]}){
	                        if($2 eq "+" && $tran_info{$arr[6]} eq "+"){
	                                $intron = "intron$exon";
	                        }elsif($2 eq "-" && $tran_info{$arr[6]} eq "+"){
        	                        $exon--;
	                                $intron = "intron$exon";
	                        }elsif($2 eq "+" && $tran_info{$arr[6]} eq "-"){
        	                        $exon--;
                	                $intron = "intron$exon";
	                        }elsif($2 eq "-" && $tran_info{$arr[6]} eq "-"){
        	                        $exon--;$exon--;
	                                $intron = "intron$exon";
        	                }
	                }
	                $annotation = "$arr[6]\t$tran\t$intron\t$nt_change\tnil";
	        }elsif($arr[5] eq "exonic" && $arr[8] ne "."){
			my @crr = split/,/,$arr[9];
        		my ($gene,$tran,$region,$nt_change,$aa_change) = split/:/,$crr[0];
		        for my $anno(@crr){
        		        my @drr = split/:/,$anno;
                		next unless(exists $genetran{$drr[0]} && $genetran{$drr[0]} eq $drr[1]);
		                ($gene,$tran,$region,$nt_change,$aa_change) = split/:/,$anno;
		        }
	        	$annotation = "$gene\t$tran\t$region\t$nt_change\t$aa_change";
		}elsif($arr[5] eq "UTR3" && $arr[7] ne "."){
	                my @crr = split/,/,$arr[7];
	                my ($tran,$nt_change) = split/:/,$crr[0];
	                for my $anno(@crr){
	                        my @drr = split/:/,$anno;
	                        next unless(exists $genetran{$arr[6]} && $genetran{$arr[6]} eq $drr[0]);
	                        ($tran,$nt_change) = split/:/,$anno;
	                }
	                $annotation = "$arr[6]\t$tran\t.\t$nt_change\tnil";
	        }else{
	                $annotation = "$arr[6]\t.\t.\t.\tnil";
	        }
	}
	print FLS "$Mut\t$dbSNP\t$annotation\t$SIFT\t$clinvar\t$cosmid\t$cosmic_occu\t$EPS6500_AF\t$kgenome_AF\t$tumor_AF\t$tumor_altread\t$Multiple_Overlap_ref\t$Multiple_Single_ref\t$One_Overlap_ref\t$One_Single_ref\t$Multiple_Overlap_alt\t$Multiple_Single_alt\t$One_Overlap_alt\t$One_Single_alt\t$tumor_depth\t$COUNT\t$MEAN_SUPPORT\t$MEAN_FREQUENCY\n";
}
close FL;
close FLS;
