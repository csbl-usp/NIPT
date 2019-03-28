#!/usr/bin/perl -w

#Author: Jaqueline Wang
#MsC in Bioinformatics Graduate Program - USP

#SCRIPT TO SEPARATE READ WITH GOOD MAPPING

#The program receives 4 parameters:
#1 - Bam file
#2 - Mapping quality
#3 - Trio number
#4 - Sample type

####################################################################################

use strict;
use Getopt::Long;

my $help = 0;
my $BAM;
my $qual = 20;
my $amostra;
my $trio;

GetOptions("help|h" => \$help,
	   "b=s" => \$BAM,
	   "m=s" => \$qual,
	   "t=s" => \$trio,
	   "a=s" => \$amostra
    ) or die "Failed to take the options! \n";

if ($help || !($BAM && $trio && $amostra)) {die "\
TThis script requires three inputs, the bam file, the trio number and the sample type. \
Other parameters have default values, but can be changed. \
The outputs are two files. \
	BOM_MAP - contain the reads that have a mapping quality higher than the threshold. \ 
	MAU_MAP - contain the reads that have a mapping quality lower than the threshold. \
\
Required parameters: \
	-b	Bam files \
	-t	Trio number \
	-a	Type of sample AF (alleged father), M (mother) ou P (plasma) \ 
\
Other parameters: \
	-h	Show the options \
	-m	Mapping quality of reads (default = 20) \
\n";
}

####################################################################################

my $nome = "Trio$trio"."_Sample$amostra";

#Armazenamos os cromossomos e intervalos escolhidos como micro-haplótipos
open (POS, "Files/Microhaplotypes.txt") or die "Failed to obtain the microhaplotypes! \n";

while (my $pos = <POS>) {
    chomp ($pos);

    #Armazenamos o cromossomo e a posição inicial e final do microhaplótipo
    $pos =~ m/(chr[0-9]+):(.*)-(.*)/;
    my $chr = $1;
    my $ini = $2;
    my $fim = $3;

    #Abrimos o arquivo BOM_ALN e extraímos as informações de mapeamento
    open (INFILE, "$nome.BOM_ALN/$nome.$chr.$ini-$fim.BOM_ALN.tsv") or die "Failed to open BOM_ALN file! \n";

    #Mapeamento BOM
    open (OUT1, ">$nome.$chr.$ini-$fim.BOM_MAP.M$qual.tsv") or die "Failed to open BOM_MAP file! \n";

    #Mapeamento RUIM
    open (OUT2, ">$nome.$chr.$ini-$fim.MAU_MAP.M$qual.tsv") or die "Failed to open MAU_MAP file! \n";
    
    while (my $line1 = <INFILE>) {
	my @array = split (/\t/, $line1);

	if ($array[4] >= $qual) {
	    print OUT1 $line1;
	}

	else {
	    print OUT2 $line1;
	}
    } # while (my $line1 = <INFILE>)

    close (INFILE);
    close (OUT1);
    close (OUT2);
} #while (my $pos = <POS>)

close (POS);
system ("mkdir $nome.MAU_MAP.M$qual $nome.BOM_MAP.M$qual");
system ("mv $nome.*.MAU_MAP.M$qual.tsv $nome.MAU_MAP.M$qual");
system ("mv $nome.*.BOM_MAP.M$qual.tsv $nome.BOM_MAP.M$qual");
