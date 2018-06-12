# This program is for single tumor ffpe sample data analysis

GMP_ConSeq_se.pl is main program

1. Fastq clean raw data
2. BWA map
3. Mutscan scan hotspot mutation
4. Genefuse find fusion
5. samtools sort 
6. ConSeqV2.4 or Duplex_sam remove PCR dup and build consensus sequence
7. Sort sam and mpileup file
8. VarScan2 call SNV/Indel
9. Annovar annotation
10. MrBam check mutation read map
11. Baseline annotation
12. Extract tab file 
13. Four Hotspot filter(RGFR L858R、T790M、19del、E13-A20 fusion)
