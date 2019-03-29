#!/usr/bin/perl -w

#Author: Jaqueline Wang
#MsC in Bioinformatics Graduate Program - USP

#SCRIPT TO CALCULATE THE PROBABILITY OF PATERNITY USING BAM FILES FROM ALLEGED FATHER, MOTHER AND PLASMA.

#The program requires 3 bam files
#S = Alleged father
#M = Mother
#P = Plasma

#The program requires 8 parameters
#1 - Mapping quality
#2 - YesCIGAR ou NotCIGAR
#3 - Bases quality
#4 - Percentage of covered bases
#5 - Coverage of regions in alleged father and mother
#6 - Coverage of regions in plasma
#7 - 1000G population
#8 - Trio number

#Inbalance parameters
#1 - Superior
#2 - Errors
#3 - Doubt between three

####################################################################################

use strict;
use Getopt::Long;

my $help = 0;
my $BAM_SP;
my $BAM_P;
my $BAM_M;
my $cobPL = 1000;
my $map = 20;
my $cigar = "NotCIGAR";
my $qual = 20;
my $por = 70;
my $cob = 20;
my $erros = 10;
my $superior = 80;
my $duvida = 35;
my $Genome = "ALL";
my $trio;

GetOptions("help|h" => \$help,
	   "F=s" => \$BAM_SP,
	   "M=s" => \$BAM_M,
	   "P=s" => \$BAM_P,
	   "m=s" => \$map,
	   "l=s" => \$cigar,
	   "q=s" => \$qual,
	   "p=s" => \$por,
	   "c=s" => \$cob,
	   "f=s" => \$cobPL,
	   "e=s" => \$erros,
	   "s=s" => \$superior,
	   "d=s" => \$duvida,
	   "g=s" => \$Genome,
	   "t=s" => \$trio
    ) or die "Failed to take the options! \n";

if ($help || !($BAM_SP && $BAM_M && $BAM_P && $trio)) {die "\
This script is used to calculate the probability of paternity. \
It requires four inputs, the alleged father bam file, the mother bam file, the plasma bam file and the trio number. \
Other parameters have default values, but can be changed. \
\
Required parameters: \
	-F	ALLEGED FATHER bam file \
	-M	MOTHER bam file \
	-P	PLASMA bam file \
	-t	Trio number \
\
Other parameters: \
	-h	Show the options \
	-m	Mapping quality of reads (default = 20) \
	-l	Use or not the CIGAR info. If NOT (NotCIGAR), consider only (mis)matches. Other option is YesCIGAR (default = NotCIGAR) \
	-q	Bases quality (default = 20) \
	-p	Percentage of covered bases (default = 70) \
	-c	Coverage of ALLEGED FATHER and MOTHER (default = 20) \
	-f	Coverage of PLASMA (default = 1000) \
	-g	Population of 1000 Genomes (default = ALL) \
	-e	Limit for sequencing errors (default = 10) \
	-s	Limit to consider HOMOZYGOUS (default = 80) \
	-d	Limit to consider HETEROZYGOUS when ther are more than 3 possibilities (default = 35) \
\n";
}

####################################################################################

system("mkdir BAM/Trio$trio");

#Analysis of ALLEGED FATHER
system("Scripts_PaternityCalc/a_analyze_sample.pl -t $trio -b $BAM_SP -m $map -l $cigar -q $qual -p $por -c $cob -a AF -e $erros -s $superior -d $duvida");

#Analysis of MOTHER
system("Scripts_PaternityCalc/a_analyze_sample.pl -t $trio -b $BAM_M -m $map -l $cigar -q $qual -p $por -c $cob -a M -e $erros -s $superior -d $duvida");

#Analysis of PLASMA
system("Scripts_PaternityCalc/a_analyze_sample.pl -t $trio -b $BAM_P -m $map -l $cigar -q $qual -p $por -c $cob -a P -e $erros -s $superior -d $duvida");

#PROBABILITY OF PATERNITY calculation
system("Scripts_PaternityCalc/b_calculate_W.pl -t $trio -b $BAM_P -m $map -l $cigar -q $qual -p $por -c $cob -f $cobPL -g $Genome -e $erros -s $superior -d $duvida > Paternity_Calc.M$map.$cigar.Q$qual.P$por.C$cob.CP$cobPL.$Genome.E$erros.S$superior.D$duvida.txt");
