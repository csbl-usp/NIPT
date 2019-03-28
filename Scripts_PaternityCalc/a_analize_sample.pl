#!/usr/bin/perl -w

#Author: Jaqueline Wang
#MsC in Bioinformatics Graduate Program - USP

#SCRIPT TO EXTRACT THE ALLEGED FATHER AND MOTHER GENOTYPES AND OBTAIN THE PLASMA HAPLOTYPES.

#It uses the following scripts:
#1_filter_alignment.pl                         
#2_filter_mapping.pl               
#3_filter_coverage_microhaplotype_yes_cigar.pl
#3_filter_coberage_microhaplotype_not_cigar.pl
#4_filter_base_quality.pl           
#5_haplotype_freq_micro.pl             
#6_extract_genotype.pl
#7_make_report.pl


#The program requires 7 parameters:
#1 - Bam file
#2 - Mapping quality
#3 - YesCIGAR ou NotCIGAR
#4 - Bases quality
#5 - Percentage of covered bases
#6 - Coverage of regions
#7 - Sample type

#Inbalance parameters
#1 - Superior
#2 - Errors
#3 - Doubt between three

####################################################################################

use strict;
use Getopt::Long;

my $help = 0;
my $BAM;
my $map = 20;
my $cigar = "NotCIGAR";
my $qual = 20;
my $por = 70;
my $cob = 20;
my $erros = 10;
my $superior = 80;
my $duvida = 35;
my $amostra;

GetOptions("help|h" => \$help,
	   "b=s" => \$BAM,
	   "m=s" => \$map,
	   "l=s" => \$cigar,
	   "q=s" => \$qual,
	   "p=s" => \$por,
	   "c=s" => \$cob,
	   "e=s" => \$erros,
	   "s=s" => \$superior,
	   "d=s" => \$duvida,
	   "a=s" => \$amostra
    ) or die "Failed to take the options! \n";

if ($help || !($BAM && $amostra)) {die "\
This script receives the parameters to analyze the sample data (alleged father, mother or plasma). \
The output are the sample genotypes and a quality report for the data. \
\
Parameters: \
	-h	Show the options \
	-b	BAM file \
	-m	Mapping quality of reads (default = 20) \
	-l	YesCIGAR or NotCIGAR (default = NotCIGAR) \
	-q	Bases quality (default = 20) \
	-p	Percentage of covered bases (default = 70) \
	-c	Coveraged reads (default = 20) \
	-a	Type of sample AF (alleged father), M (mother) ou P (plasma) \
	-e	Limit for sequencing errors (default = 10) \
	-s	Limit to consider HOMOZYGOUS (default = 80) \
	-d	Limit to consider HETEROZYGOUS when there are more than 3 possibilities (default = 35) \
\n";
}

my $nome_genotipo = "Genotypes_$amostra.M$map.$cigar.Q$qual.P$por.C$cob.E$erros.S$superior.D$duvida.txt";
my $nome_report = "MCoverage_$amostra.M$map.$cigar.Q$qual.P$por.C$cob.E$erros.S$superior.D$duvida.txt";

#Usamos esse match para armazenar o início do nome de cada arquivo
$BAM =~ m/(IonXpress[\._]{1}([0-9]{3}))[\._]{1}R[\._]{1}([0-9]{4}_[0-9]{2}_[0-9]{2})((.)*).bam/;
my $nome = $1;
my $id = $2;
my $data = $3;

##Roda todos os scripts para análise

system ("Scripts_PaternityCalc/1_filter_alignment.pl -b $BAM");

system ("Scripts_PaternityCalc/2_filter_mapping.pl -b $BAM -m $map");

if ($cigar eq "NotCIGAR") {
    system ("Scripts_PaternityCalc/3_filter_coverage_microhaplotype_not_cigar.pl -b $BAM -m $map");
}

elsif ($cigar eq "YesCIGAR") {
    system ("Scripts_PaternityCalc/3_filter_coverage_microhaplotype_yes_cigar.pl -b $BAM -m $map");
}
    
system ("Scripts_PaternityCalc/4_filter_base_quality.pl -b $BAM -m $map -l $cigar -q $qual -p $por");
   
system ("Scripts_PaternityCalc/5_haplotype_freq_micro.pl -b $BAM -m $map -l $cigar -q $qual -p $por");

system ("mkdir $nome");

system ("mv $nome.* $nome");

if (($amostra eq "AF") || ($amostra eq "M")) {
    system ("Scripts_PaternityCalc/6_extract_genotype.pl -b $BAM -m $map -l $cigar -q $qual -p $por -c $cob -e $erros -s $superior -d $duvida > $nome_genotipo");
}

system("Scripts_PaternityCalc/7_make_report.pl -b $BAM -m $map -l $cigar -q $qual -p $por -c $cob -a $amostra -e $erros -s $superior -d $duvida > $nome_report");
