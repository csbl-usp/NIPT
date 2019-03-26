#!/usr/bin/perl -w

#Autora: Jaqueline Wang
#Mestre do Programa Interunidades em Bioinformática - USP

#SCRIPT PARA CRIAR O META FILE DAS FREQUÊNCIAS DOS HAPLÓTIPOS DO 1000 GENOMES.

#Devem ser passados 3 parâmetros
#1 - Arquivo VCF
#2 - Lista dos SNPs que compõem o micro-haplótipo
#3 - Número do Micro-haplótipo

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

system ("/home/jaque/Desktop/Scripts_Martin/Script_1000G/a_calcula_freq_TODOS_POP_SUPER.pl -v $vcf -m $micro -n $num");

system ("/home/jaque/Desktop/Scripts_Martin/Script_1000G/b_faz_META_FILE.pl -n $num");
