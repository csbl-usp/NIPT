#!/usr/bin/perl -w

#Autora: Jaqueline Wang
#Mestre do Programa Interunidades em Bioinformática - USP

#SCRIPT PARA CRIAR OS ARQUIVOS DE SAÍDA PARA TODAS AS POPULAÇÕES E SUPER-POPULAÇÕES

#Devem ser passados 3 parâmetros
#1 - Arquivo VCF
#2 - SNPs do Micro-haplótipo
#3 - Número do Micro-haplótipo

#Endereços importantes (MARTIN)
#Lista de populações: /home/jaque/Desktop/Scripts_Martin/Arquivos/lista_pop.txt
#Script utilizado: /home/jaque/Desktop/Scripts_Martin/Script_1000G/1_calcula_frequencia.pl
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
) or die "Erro ao pegar as opções! \n";

if ($help || !($vcf && $micro && $num)) {die "\
$0: \
    Esse script recebe três arquivos, um é a lista de SNPs que compõem o microhaplótipo, o segundo é o vcf do 1000 Genomes e o terceiro é o número do micro-haplótipo.\
\
Parâmetros:\
     -h ou --help : Mostra as opções\
     -v : Arquivo vcf do 1000 Genomes\
     -m : Lista com os SNPs que compõem o microhaplótipo\
     -n : Número do micro-haplótipo\

\n";
}

####################################################################################

open (POP, "/home/jaque/Desktop/Scripts_Martin/Arquivos/lista_pop.txt") or die "Não foi possível abrir o arquivo de lista de populacoes! \n";

system ("mkdir Freq_haplotipo/M$num");

while (my $line1 = <POP>) {
    chomp ($line1);
    my ($sigla, $pop) = split(/\t/, $line1);

    if ($sigla eq "TODOS") {

	system ("/home/jaque/Desktop/Scripts_Martin/Script_1000G/1_calcula_frequencia.pl -v $vcf -m $micro > Freq_haplotipo/M$num/M$num.$pop.freq.txt");
    }

    elsif ($sigla eq "s") {

	system ("/home/jaque/Desktop/Scripts_Martin/Script_1000G/1_calcula_frequencia.pl -v $vcf -m $micro -s $pop > Freq_haplotipo/M$num/M$num.$pop.freq.txt");
    }

    elsif ($sigla eq "p") {

	system ("/home/jaque/Desktop/Scripts_Martin/Script_1000G/1_calcula_frequencia.pl -v $vcf -m $micro -p $pop > Freq_haplotipo/M$num/M$num.$pop.freq.txt");
    }


}

close (POP);
