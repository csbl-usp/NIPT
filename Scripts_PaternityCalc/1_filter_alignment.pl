#!/usr/bin/perl -w

#Author: Jaqueline Wang
#MsC in Bioinformatics Graduate Program - USP

#SCRIPT TO SEPARATE READ WITH A GOOD ALIGNMENT 

#The program receives 3 parmeters:
#1 - Bam file
#2 - Trio Number
#3 - Sample type

####################################################################################

use strict;
use Getopt::Long;

my $help = 0;
my $BAM;
my $trio;
my $amostra;

GetOptions("help|h" => \$help,
	   "b=s" => \$BAM,
	   "t=s" => \$trio,
	   "a=s" => \$amostra
    ) or die "Failed to take the options! \n";

if ($help || !($BAM && $trio && $amostra)) {die "\
This script receives the bam file, the trio number and the sample type. \
The outputs are two files with sequences in bam format. \
	BOM_ALN - contain reads aligned only in ONE region. \
	MAU_ALN - contain reads aligned in more than ONE region. \
\
Parameters: \
	-h	Show the options \
	-b	Bam file \
	-a	Type of sample AF (alleged father), M (mother) ou P (plasma) \ 
	-t	Trio number \
\n";
}

####################################################################################

my $nome = "Trio$trio"."_Sample$amostra";

#Armazenamos os cromossomos e intervalos escolhidos como micro-haplótipos
 
open (POS, "Files/Microhaplotypes.txt") or die "Failed to obtain the microhaplotypes! \n";

#Indexamos o BAM para poder visualizar os dados
system ("Samtools/samtools-1.3.1/samtools index $BAM");

while (my $pos = <POS>){

    chomp ($pos);

    #Armazenamos o cromossomo e a posição inicial e final do microhaplótipo
    $pos =~ m/(chr[0-9]+):(.*)-(.*)/;
    my $chr = $1;
    my $ini = $2;
    my $fim = $3;

    #Extraímos as informações do BAM e colocar em um arquivo TSV
    system ("Samtools/samtools-1.3.1/samtools view $BAM $pos > $nome.$chr.$ini-$fim.tsv");

    #Abrimos o aquivo TSV e extraímos as informações das regiões analisadas
    open (TSV, "$nome.$chr.$ini-$fim.tsv") or die "Failed to open TSV file! \n";

    #Alinhamentos em uma região
    open (OUT1, ">$nome.$chr.$ini-$fim.BOM_ALN.tsv") or die "Failed to open OUT1 file! \n";

    #Alinhamentos em duas ou mais regiões
    open (OUT2, ">$nome.$chr.$ini-$fim.MAU_ALN.tsv") or die "Failed to open OUT2 file! \n";

    while (my $line = <TSV>){
	my @array = split (/\t/, $line);
    
	if ((($array[1] == 16) || ($array[1] == 0))) {
	    print OUT1 $line;
	} 

	else {
	    print OUT2 $line;
	}
    } #while (my $line = <TSV>)

    close (TSV);
    close (OUT1);
    close (OUT2);
    system ("rm $nome.$chr.$ini-$fim.tsv");
} #while (my $pos = <POS>)

close (POS);
system ("rm $nome*.bai");
system ("mkdir $nome.MAU_ALN $nome.BOM_ALN");
system ("mv $nome.*.MAU_ALN.tsv $nome.MAU_ALN");
system ("mv $nome.*.BOM_ALN.tsv $nome.BOM_ALN");
