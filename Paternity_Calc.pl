#!/usr/bin/perl -w

#Autora: Jaqueline Wang
#Mestre do Programa Interunidades em Bioinformática - USP

#SCRIPT PARA CALCULAR A PROBABILIDADE DE PATERNIDADE UTILIZANDO OS ARQUIVOS DO SUPOSTO PAI, DA MÃE E DO PLASMA.

#Devem ser passados 3 arquivos BAM
#S = Suposto Pai
#M = Mãe
#P = Plasma

#Devem ser passados 7 parâmetros:
#1 - Qualidade do mapeamento
#2 - ComCIGAR ou SemCIGAR
#3 - Qualidade das BASES
#4 - Porcentagem das BASES cobertas
#5 - Cobertura das regiões do SUPOSTO PAI e da MÃE
#6 - Cobertura das regiões do PLASMA
#7 - População do Banco de Dados

#Parâmetros de desbalanço
#1 - Superior
#2 - Erros
#3 - Dúvida de três

#Endereços importantes (MARTIN)
#Micro-haplótipos: /home/jaque/Desktop/Scripts_Martin/Arquivos/microhaplotipos.txt

####################################################################################

use strict;
use Getopt::Long;

my $help = 0;
my $BAM_SP;
my $BAM_P;
my $BAM_M;
my $cobPL = 1000;
my $map = 20;
my $cigar = "SemCIGAR";
my $qual = 20;
my $por = 70;
my $cob = 20;
my $erros = 10;
my $superior = 80;
my $duvida = 35;
my $Genome = "TODOS";

GetOptions("help|h" => \$help,
	   "X=s" => \$BAM_SP,
	   "Y=s" => \$BAM_M,
	   "Z=s" => \$BAM_P,
	   "m=s" => \$map,
	   "l=s" => \$cigar,
	   "q=s" => \$qual,
	   "p=s" => \$por,
	   "c=s" => \$cob,
	   "f=s" => \$cobPL,
	   "e=s" => \$erros,
	   "s=s" => \$superior,
	   "d=s" => \$duvida,
	   "g=s" => \$Genome
    ) or die "Erro ao pegar as opções! \n";

if ($help || !($BAM_SP && $BAM_M && $BAM_P)) {die "\
$0: \
Esse script recebe os parâmetros para análise dos dados e a cobertura mínima do plasma. A saída são os genótipos encontrados para a amostra do SUPOSTO PAI, da MÃE e do PLASMA. Para cada micro-haplótipo é obtida a Evidência da Paternidade, e posteriormente calculado o Índice de Paternidade (IP) e a Probabilidade de Paternidade Posterior (W).\
\
Parâmetros:\
    -h ou --help : Mostra as opções\
    -S : Arquivo BAM do SUPOSTO PAI\
    -M : Arquivo BAM da MÃE\
    -P : Arquivo BAM do PLASMA\
    -m : Qualidade do mapeamento dos reads (padrão = 20)\
    -l : ComCIGAR ou SemCIGAR (padrão = SemCIGAR)
    -q : Qualidade das bases (padrão = 20)\
    -p : Porcentagem de bases cobertas (padrão = 70)\
    -c : Cobertura de reads do SUPOSTO PAI e MÃE (padrão = 20)\
    -f : Cobertura de reads do PLASMA (padrão = 1000)\
    -g : Populacão do 1000 Genomes (padrão = TODOS)\
    -e : Limite para erros de sequenciamento (Padrão = 10)\
    -s : Limite para considerar HOMOZIGOTO (Padrão = 80)\
    -d : Limite para considerar HETEROZIGOTO quando existem 3 ou mais (Padrão = 35)\
\n";
}

####################################################################################

#Fazemos a análise do arquivo do SUPOSTO PAI
system("/home/jaque/Desktop/Scripts_Martin/a_analisa_amostra.pl -b $BAM_SP -m $map -l $cigar -q $qual -p $por -c $cob -a SP -e $erros -s $superior -d $duvida");

#Fazemos a análise do arquivo da MÃE
system("/home/jaque/Desktop/Scripts_Martin/a_analisa_amostra.pl -b $BAM_M -m $map -l $cigar -q $qual -p $por -c $cob -a M -e $erros -s $superior -d $duvida");

#Fazemos a análise do arquivo do PLASMA
system("/home/jaque/Desktop/Scripts_Martin/a_analisa_amostra.pl -b $BAM_P -m $map -l $cigar -q $qual -p $por -c $cob -a P -e $erros -s $superior -d $duvida");

#Fazemos o cálculo da PROBABILIDADE DE PATERNIDADE
system("/home/jaque/Desktop/Scripts_Martin/b_calcula_W.pl -b $BAM_P -m $map -l $cigar -q $qual -p $por -c $cob -f $cobPL -g $Genome -e $erros -s $superior -d $duvida > Paternity_Calc.M$map.$cigar.Q$qual.P$por.C$cob.CP$cobPL.$Genome.E$erros.S$superior.D$duvida.txt");

    
    
