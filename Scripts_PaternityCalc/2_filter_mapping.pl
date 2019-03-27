#!/usr/bin/perl -w

#Author: Jaqueline Wang
#MsC in Bioinformatics Graduate Program - USP

#SCRIPT TO SEPARATE READ WITH GOOD MAPPING

#The program receives 2 parameters:
#1 - Bam file
#2 - Mapping quality

####################################################################################

use strict;
use Getopt::Long;

my $help = 0;
my $BAM;
my $qual = 20;

GetOptions("help|h" => \$help,
	   "b=s" => \$BAM,
	   "m=s" => \$qual
    ) or die "Erro ao pegar as opções! \n";

if ($help || !($BAM)) {die "\
This script receives two inputs, the bam file and the mapping quality. \
The outputs are two files with the sequences in bam format. \
BOM_MAP - contain the reads that have a mapping quality higher than the threshold. \ 
MAU_MAP - contain the reads that have a mapping quality lower than the threshold. \
\
Parameters:\
    -h ou --help : Show the options\
    -b : Bam files\
    -m : Mapping quality of reads (default = 20)\
\n";
}

####################################################################################

#Usamos esse match para armazenar o início do nome de cada arquivo
$BAM =~ m/(IonXpress[\._]{1}([0-9]{3}))[\._]{1}R[\._]{1}([0-9]{4}_[0-9]{2}_[0-9]{2})((.)*).bam/;

my $nome = $1;
my $id = $2;
my $data = $3;

#Armazenamos os cromossomos e intervalos escolhidos como micro-haplótipos
open (POS, "Files/microhaplotypes.txt") or die "Failed to obtain the microhaplotypes! \n";

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
