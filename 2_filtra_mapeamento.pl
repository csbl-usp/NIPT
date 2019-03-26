#!/usr/bin/perl -w

#Autora : Jaqueline Wang
#Mestre do programa Interunidades em Bioinformática - USP

#SCRIPT PARA SEPARAR READS QUE TENHAM UM BOM MAPEAMENTO.

#Devem ser passados 2 parâmetros:
#1 - Arquivo BAM
#2 - Qualidade do mapeamento

#Endereços importantes (MARTIN)
#Micro-haplótipos: /home/jaque/Desktop/Scripts_Martin/Arquivos/microhaplotipos.txt

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
$0: \
Esse script recebe duas entradas, o arquivo BAM e a qualidade de mapeamento. A saída são dois arquivos com as sequências no formato BAM: \
BOM_MAP - contém os reads que estejam com um mapeamento com qualidade maior que o desejado. \ 
MAU_MAP - contém os reads que estejam com um mapeamento com qualidade menor que o desejado. \
\
Parâmetros:\
    -h ou --help : Mostra as opções\
    -b : Arquivo BAM\
    -m : Qualidade do mapeamento dos reads (padrão = 20)\
\n";
}

####################################################################################

#Usamos esse match para armazenar o início do nome de cada arquivo
$BAM =~ m/(IonXpress[\._]{1}([0-9]{3}))[\._]{1}R[\._]{1}([0-9]{4}_[0-9]{2}_[0-9]{2})((.)*).bam/;

my $nome = $1;
my $id = $2;
my $data = $3;

#Armazenamos os cromossomos e intervalos escolhidos como micro-haplótipos
open (POS, " /home/jaque/Desktop/Scripts_Martin/Arquivos/microhaplotipos.txt") or die "Não foi possível obter os micro-haplótipos! \n";

while (my $pos = <POS>) {
    
    chomp ($pos);

    #Armazenamos o cromossomo e a posição inicial e final do microhaplótipo
    $pos =~ m/(chr[0-9]+):(.*)-(.*)/;
    my $chr = $1;
    my $ini = $2;
    my $fim = $3;

    #Abrimos o arquivo BOM_ALN e extraímos as informações de mapeamento
    open (INFILE, "$nome.BOM_ALN/$nome.$chr.$ini-$fim.BOM_ALN.tsv") or die "Não foi possível abrir o arquivo com os BOMs ALINHAMENTOS! \n";

    #Mapeamento BOM
    open (OUT1, ">$nome.$chr.$ini-$fim.BOM_MAP.M$qual.tsv") or die "Não foi possível abrir o arquivo com BOM MAPEAMENTO! \n";

    #Mapeamento RUIM
    open (OUT2, ">$nome.$chr.$ini-$fim.MAU_MAP.M$qual.tsv") or die "Não foi possível abrir o arquivo com MAU MAPEAMENTO! \n";
    
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
