#!/usr/bin/perl -w

#Author: Jaqueline Wang
#MsC in Bioinformatics Graduate Program - USP

#SCRIPT TO CREATE A META FILE OF HAPLOTYPES FREQUENCIES FROM 1000 GENOMES.

#The script requires 3 parameters
#1 - Vcf file
#2 - List of SNPs that compose the microhaplotype
#3 - Microhaplotype number

####################################################################################

use strict;
use Getopt::Long;

my $help = 0;
my ($vcf, $micro);
my $num = 0;

GetOptions("help|h" => \$help,
	   "v=s" => \$vcf,
	   "m=s" => \$micro,
	   "n=s" => \$num,
) or die "Failed to take the options! \n";

if ($help || !($vcf && $micro && $num)) {die "\
This script requires the following three parameters. \
\
Parameters: \
     -h	Show the options \
     -v	Vcf file from 1000 Genomes \
     -m	List of SNPs that compose the microhaplotype \
     -n	Microhaplotype number \
\n";
}

####################################################################################

system ("Scripts_1000G/a_calc_freq_ALL_POP_SUPER.pl -v $vcf -m $micro -n $num");

system ("Scripts_1000G/b_make_META_FILE.pl -n $num");
