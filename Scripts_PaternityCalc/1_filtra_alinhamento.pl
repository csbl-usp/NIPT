#!/usr/bin/perl -w

#Autora: Jaqueline Wang
#Mestre do Programa Interunidades em Bioinformática - USP

#SCRIPT PARA SEPARAR READS QUE TENHAM UM BOM ALINHAMENTO.

#Deve ser passado 1 parâmetro:
#1 - Arquivo BAM

#Endereços importantes (MARTIN)
#Samtools: /home/jaque/Desktop/Scripts_Martin/Samtools/samtools-1.3.1/samtools
#Micro-haplótipos: /home/jaque/Desktop/Scripts_Martin/Arquivos/microhaplotipos.txt

####################################################################################

use strict;
use Getopt::Long;

my $help = 0;
my $BAM;

GetOptions("help|h" => \$help,
	   "b=s" => \$BAM
    ) or die "Erro ao pegar as opções! \n";

if ($help || !($BAM)) {die "\
$0: \
Esse script recebe uma entrada, o arquivo BAM. A saída são dois arquivos com as sequências no formato BAM: \
BOM_ALN - contém reads que estejam alinhados em uma região apenas. \
MAU_ALN - contém reads que estejam alinhados ema mais de uma região. \
\
Parâmetros:\
    -h ou --help : Mostra as opções\
    -b : Arquivo BAM\
\n";
}

####################################################################################

#Usamos esse match para armazenar o início do nome de cada arquivo
$BAM =~ m/(IonXpress[\._]{1}([0-9]{3}))[\._]{1}R[\._]{1}([0-9]{4}_[0-9]{2}_[0-9]{2})((.)*).bam/;

my $nome = $1;
my $id = $2;
my $data = $3;

#Armazenamos os cromossomos e intervalos escolhidos como micro-haplótipos
 
open (POS, "/home/jaque/Desktop/Scripts_Martin/Arquivos/microhaplotipos.txt") or die "Não foi possível obter os micro-haplótipos! \n";

#Indexamos o BAM para poder visualizar os dados
system ("/home/jaque/Desktop/Scripts_Martin/Samtools/samtools-1.3.1/samtools index $BAM");

while (my $pos = <POS>){

    chomp ($pos);

    #Armazenamos o cromossomo e a posição inicial e final do microhaplótipo
    $pos =~ m/(chr[0-9]+):(.*)-(.*)/;
    my $chr = $1;
    my $ini = $2;
    my $fim = $3;

    #Extraímos as informações do BAM e colocar em um arquivo TSV
    system ("/home/jaque/Desktop/Scripts_Martin/Samtools/samtools-1.3.1/samtools view $BAM $pos > $nome.$chr.$ini-$fim.tsv");

    #Abrimos o aquivo TSV e extraímos as informações das regiões analisadas
    open (TSV, "$nome.$chr.$ini-$fim.tsv") or die "Não foi possível abrir o arquivo TSV! \n";

    #Alinhamentos em uma região
    open (OUT1, ">$nome.$chr.$ini-$fim.BOM_ALN.tsv") or die "Não foi possível abrir o arquivo OUT1! \n";

    #Alinhamentos em duas ou mais regiões
    open (OUT2, ">$nome.$chr.$ini-$fim.MAU_ALN.tsv") or die "Não foi possível abrir o arquivo OUT2! \n";

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
