#!/usr/bin/perl -w
use strict;
use File::Basename;
use FindBin '$Bin';

my ($outdir,$samplename,$cftype,$cfDNAfq1,$cfDNAfq2,$bed,$t,$p5index) = @ARGV;
=head1 Uasge

        perl  ctDNA_pairV2.pl <outdir> <samplename> <cfDNAtype> <cfDNAfq1> <cfDNAfq2> <bed> <thread> <p5index>
	
	<outdir>		output dir
	<samplename>		sample  prefix
	<cfDNAtype>		Index/Duplex/nonUMI
	<cfDNAfq1>		cfDNA fq1
	<cfDNAfq2>              cfDNA fq2
	<bed>			bed file
	<thread>		thread use
	<p5index>		Index2 list remove

=cut

die `pod2text $0` if (@ARGV != 8);

########################cfDNA#######################
system("mkdir -m 755 -p $outdir") if (!-d "$outdir");
system("mkdir -m 755 -p $outdir/mutscan") if (!-d "$outdir/mutscan");
system("mkdir -m 755 -p $outdir/fusion") if (!-d "$outdir/fusion");

open OUT1,">$outdir/$samplename\_css\_cfDNA.sh";
print OUT1 "echo;date\n";
if($cftype eq "Index"){
	print OUT1 "$Bin/fastp -f 1 -F 1 -t 1 -T 1 -3 -W 4 -M 25 --filter_by_index2=$p5index --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -i $cfDNAfq1 -o $outdir/$samplename\_cfDNA_good_R1.fastq.gz -I $cfDNAfq2 -O $outdir/$samplename\_cfDNA_good_R2.fastq.gz -U --umi_loc=index2 --umi_prefix=UMI -h $outdir/$samplename\_cfDNA_fastp.html -j $outdir/$samplename\_cfDNA_fastp.json > $outdir/$samplename\_cfDNA_fastp.stat\n";
}elsif($cftype eq "nonUMI"){
	print OUT1 "$Bin/fastp -f 1 -F 1 -t 1 -T 1 -3 -W 4 -M 25 --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -i $cfDNAfq1 -o $outdir/$samplename\_cfDNA_good_R1.fastq.gz -I $cfDNAfq2 -O $outdir/$samplename\_cfDNA_good_R2.fastq.gz  -h $outdir/$samplename\_cfDNA_fastp.html -j $outdir/$samplename\_cfDNA_fastp.json > $outdir/$samplename\_cfDNA_fastp.stat\n";
}elsif($cftype eq "Duplex"){
	print OUT1 "$Bin/fastp --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -i $cfDNAfq1 -o $outdir/$samplename\_cfDNA_good_R1.fastq.gz -I $cfDNAfq2 -O $outdir/$samplename\_cfDNA_good_R2.fastq.gz  -h $outdir/$samplename\_cfDNA_fastp.html -j $outdir/$samplename\_cfDNA_fastp.json  -U --umi_loc=per_read --umi_skip=4 --umi_prefix=UMI  --umi_len=3  > $outdir/$samplename\_cfDNA_fastp.stat\n";
}
print OUT1 "nohup $Bin/mutscan -m $Bin/GMP_alt.csv -S 1 -1 $outdir/$samplename\_cfDNA_good_R1.fastq.gz -2 $outdir/$samplename\_cfDNA_good_R2.fastq.gz -h $outdir/mutscan/$samplename\_mutscan_alt.html -s >$outdir/mutscan/$samplename\_mutscan.log 2>&1 &\n";
print OUT1 "nohup $Bin/mutscan -m $Bin/GMP_ref.csv -S 1 -1 $outdir/$samplename\_cfDNA_good_R1.fastq.gz -2 $outdir/$samplename\_cfDNA_good_R2.fastq.gz -h $outdir/mutscan/$samplename\_mutscan_ref.html -s >$outdir/mutscan/$samplename\_mutscan.log 2>&1 &\n";
print OUT1 "nohup $Bin/GeneFuse/genefuse -t $t -u 1 -r $Bin/WES_ref/hg19.fa -1 $outdir/$samplename\_cfDNA_good_R1.fastq.gz -2 $outdir/$samplename\_cfDNA_good_R2.fastq.gz -f $Bin/GeneFuse/genes/ALK.hg19.csv -h $outdir/fusion/$samplename\_fusion.html >$outdir/fusion/$samplename.fusion.log 2>&1 &\n";
print OUT1 "$Bin/bwa/bwa mem -k 32 -R \"\@RG\\tID:$samplename\\tLB:$samplename\\tPL:ILLUMINA\\tSM:$samplename\" -t $t -M $Bin/WES_ref/hg19.fa $outdir/$samplename\_cfDNA_good_R1.fastq.gz $outdir/$samplename\_cfDNA_good_R2.fastq.gz > $outdir/$samplename.cfDNA.sam\n";
print OUT1 "$Bin/samtools-1.3/bin/samtools sort -@ $t -T $outdir/$samplename.cfDNA -o $outdir/$samplename.cfDNA.sort.sam $outdir/$samplename.cfDNA.sam\n";
if($cftype eq "Index"){
	print OUT1 "$Bin/ConSeqV2.4.pl -e Index -b 8 -c 1 -u 0.8 -s UMI -r $Bin/WES_ref/hg19.fa -i $outdir/$samplename.cfDNA.sort.sam -o $outdir/$samplename\_css.cfDNA.sam -a $outdir/$samplename\_all.cfDNA.sam -l $outdir/$samplename\_css.cfDNA.log\n";
}elsif($cftype eq "Duplex"){
	print OUT1 "$Bin/ConSeqV2.4.pl -e Duplex -b 6 -c 1 -u 0.8 -s UMI -r $Bin/WES_ref/hg19.fa -i $outdir/$samplename.cfDNA.sort.sam -o $outdir/$samplename\_css.cfDNA.sam -a $outdir/$samplename\_all.cfDNA.sam -l $outdir/$samplename\_css.cfDNA.log\n";
}elsif($cftype eq "nonUMI"){
	print OUT1 "$Bin/ConSeqV2.4.pl -u 0.5 -s nonUMI -r $Bin/WES_ref/hg19.fa -i $outdir/$samplename.cfDNA.sort.sam -o $outdir/$samplename\_css.cfDNA.sam -a $outdir/$samplename\_all.cfDNA.sam -l $outdir/$samplename\_css.cfDNA.log\n";
}
print OUT1 "$Bin/samtools-1.3/bin/samtools sort -@ $t -T $outdir/$samplename\_css.cfDNA -o $outdir/$samplename\_css.cfDNA.dedup.bam $outdir/$samplename\_css.cfDNA.sam\n";
print OUT1 "$Bin/samtools-1.3/bin/samtools index $outdir/$samplename\_css.cfDNA.dedup.bam\n";
print OUT1 "$Bin/samtools-0.1.19/samtools mpileup -AB -Q 25 -d 100000 -f $Bin/WES_ref/hg19.fa -l $bed  $outdir/$samplename\_css.cfDNA.dedup.bam >$outdir/$samplename\_css.cfDNA.dedup.mpileup\n";
print OUT1 "perl $Bin/stat_umi2.pl $outdir/$samplename\_cfDNA_fastp.stat $outdir/$samplename\_css.cfDNA.log $outdir/$samplename\_css.cfDNA.dedup.mpileup CSS >$outdir/$samplename.cfDNA.stat\n";
print OUT1 "echo;date\n";
close OUT1;
######################################################################


#########################SOMATIC##################################

open OUT3,">$outdir/$samplename\_css\_somatic.sh";
print OUT3 "echo;date\n";
print OUT3 "java -jar $Bin/VarScan.v2.3.8.jar mpileup2snp $outdir/$samplename\_css.cfDNA.dedup.mpileup --min-coverage 80 --min-reads2 2 --min-avg-qual 25 --min-var-freq 0.0001 --p-value 0.99 --min-freq-for-hom 90 --strand-filter 0 --output-vcf 1 --variants 1 > $outdir/$samplename\_css.snp.vcf\n";
print OUT3 "java -jar $Bin/VarScan.v2.3.8.jar mpileup2indel $outdir/$samplename\_css.cfDNA.dedup.mpileup --min-coverage 80 --min-reads2 2 --min-avg-qual 25 --min-var-freq 0.0001 --p-value 0.99 --min-freq-for-hom 90 --strand-filter 0 --output-vcf 1 --variants 1 > $outdir/$samplename\_css.indel.vcf\n";
print OUT3 "$Bin/annovar/table_annovar.pl $outdir/$samplename\_css.snp.vcf $Bin/annovar/humandb/ -buildver hg19 -out $outdir/$samplename\_css.snp -remove -protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2014oct_all,1000g2014oct_afr,1000g2014oct_eas,1000g2014oct_eur,snp138,ljb26_all,cosmic77,clinvar_20160302 -operation g,r,r,f,f,f,f,f,f,f,f,f -nastring . -vcfinput\n";
print OUT3 "$Bin/annovar/table_annovar.pl $outdir/$samplename\_css.indel.vcf $Bin/annovar/humandb/ -buildver hg19 -out $outdir/$samplename\_css.indel -remove -protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2014oct_all,1000g2014oct_afr,1000g2014oct_eas,1000g2014oct_eur,snp138,ljb26_all,cosmic77,clinvar_20160302 -operation g,r,r,f,f,f,f,f,f,f,f,f -nastring . -vcfinput\n";
print OUT3 "$Bin/paraMrBam $t $outdir/$samplename\_css.snp.hg19_multianno.txt --cfdna $outdir/$samplename\_css.cfDNA.dedup.bam -m 2 -q 25 --fast >$outdir/$samplename\_css.snp_MrBam.txt\n";
print OUT3 "$Bin/paraMrBam $t $outdir/$samplename\_css.indel.hg19_multianno.txt --cfdna $outdir/$samplename\_css.cfDNA.dedup.bam -m 2 -q 25 --fast >$outdir/$samplename\_css.indel_MrBam.txt\n";
print OUT3 "$Bin/baselineanno -i $outdir/$samplename\_css.snp_MrBam.txt -o $outdir/$samplename\_css.snp_MrBam.baseline -t 1 -redis 192.168.1.10:6379 -anno 4\n";
print OUT3 "$Bin/baselineanno -i $outdir/$samplename\_css.indel_MrBam.txt -o $outdir/$samplename\_css.indel_MrBam.baseline -t 1 -redis 192.168.1.10:6379 -anno 4\n";
print OUT3 "$Bin/extractSE.pl $outdir/$samplename\_css.snp_MrBam.baseline $outdir/$samplename\_css.snp.txt\n";
print OUT3 "$Bin/extractSE.pl $outdir/$samplename\_css.indel_MrBam.baseline $outdir/$samplename\_css.indel.txt\n";
print OUT3 "perl $Bin/stat3.pl -s $outdir/$samplename\_css.snp.txt -i $outdir/$samplename\_css.indel.txt -r $outdir/mutscan/$samplename\_mutscan_ref.html -a $outdir/mutscan/$samplename\_mutscan_alt.html -f $outdir/fusion/$samplename.fusion.log >$outdir/$samplename\_css.result.txt\n";
print OUT3 "echo;date\n";
close OUT3;

open OUT4,">$outdir/$samplename\_css\_run.sh";
print OUT4 "nohup sh $outdir/$samplename\_css\_cfDNA.sh >$outdir/$samplename\_css\_cfDNA.sh.o 2>&1 &\n";
print OUT4 "wait\n";
print OUT4 "sh $outdir/$samplename\_css\_somatic.sh >$outdir/$samplename\_css\_somatic.sh.o 2>&1\n";
close OUT4;

print "cd $outdir\n";
print "nohup sh $outdir/$samplename\_css\_run.sh >$outdir/$samplename\_css\_run.sh.o 2>&1 &\n"
