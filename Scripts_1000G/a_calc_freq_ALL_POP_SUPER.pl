#!/usr/bin/perl -w

#Author: Jaqueline Wang
#MsC in Bioinformatics Graduate Program - USP

#SCRIPT TO CREATE FILES FOR EACH POPULATION AND SUPER-POPULATION

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
This script requires three parameters. \
\
Parameters: \
	-h	Show the options \
	-v	Vcf file from 1000 Genomes \
	-m	List of SNPs that compose the microhaplotype \
	-n	Microhaplotype number \
\n";
}

####################################################################################

open (POP, "Files/pop_list.txt") or die "Failed to open the POP_LIST file! \n";

system ("mkdir 1000G_test/Freq_haplo/M$num");

while (my $line1 = <POP>) {
    chomp ($line1);
    my ($sigla, $pop) = split(/\t/, $line1);

    if ($sigla eq "ALL") {
	system ("Scripts_1000G/1_calc_freq.pl -v $vcf -m $micro > 1000G_test/Freq_haplo/M$num/M$num.$pop.freq.txt");
    }

    elsif ($sigla eq "s") {
	system ("Scripts_1000G/1_calc_freq.pl -v $vcf -m $micro -s $pop > 1000G_test/Freq_haplo/M$num/M$num.$pop.freq.txt");
    }

    elsif ($sigla eq "p") {
	system ("Scripts_1000G/1_calc_freq.pl -v $vcf -m $micro -p $pop > 1000G_test/Freq_haplo/M$num/M$num.$pop.freq.txt");
    }
}

close (POP);
