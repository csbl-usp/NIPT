#!/usr/bin/perl -w

#Autora: Jaqueline Wang
#Mestre do Programa Interunidades em Bioinformática - USP

#SCRIPT PARA EXTRAIR OS GENÓTIPOS DO SUPOSTO PAI E MAE E OBTER OS HAPLÓTIPOS DO PLASMA.
#São utilizados os seguintes scripts:
#1_filtra_alinhamento.pl                         
#2_filtra_mapeamento.pl               
#3_filtra_cobertura_microhaplotipo_com_cigar.pl
#3_filtra_cobertura_microhaplotipo_sem_cigar.pl
#4_filtra_qualidade_bases.pl           
#5_haplotipo_freq_micro.pl             
#6_extrai_genotipo.pl
#7_faz_report.pl


#Devem ser passado 7 parâmetros:
#1 - Arquivo BAM
#2 - Qualidade do mapeamento
#3 - ComCIGAR ou SemCIGAR
#4 - Qualidade das BASES
#5 - Porcentagem das BASES cobertas
#6 - Cobertura das regiões
#7 - Tipo de amostra

#Parâmetros de desbalanço
#1 - Superior
#2 - Erros
#3 - Dúvida de três

#Endereços importantes (MARTIN)
#Scripts de análise: /home/jaque/Desktop/Scripts_Martin/

####################################################################################

use strict;
use Getopt::Long;

my $help = 0;
my $BAM;
my $map = 20;
my $cigar = "SemCIGAR";
my $qual = 20;
my $por = 70;
my $cob = 20;
my $erros = 10;
my $superior = 80;
my $duvida = 35;
my $amostra;

GetOptions("help|h" => \$help,
	   "b=s" => \$BAM,
	   "m=s" => \$map,
	   "l=s" => \$cigar,
	   "q=s" => \$qual,
	   "p=s" => \$por,
	   "c=s" => \$cob,
	   "e=s" => \$erros,
	   "s=s" => \$superior,
	   "d=s" => \$duvida,
	   "a=s" => \$amostra
    ) or die "Erro ao pegar as opções! \n";

if ($help || !($BAM && $amostra)) {die "\
$0: \
Esse script recebe os parâmetros para análise dos dados da amostra (SUPOSTO PAI, MÃE ou PLASMA). A saída são os genótipos encontrados para a amostra do analisada e um report da qualidade do arquivo BAM.\
\
Parâmetros:\
    -h ou --help : Mostra as opções\
    -b : Arquivo BAM\
    -m : Qualidade do mapeamento dos reads (padrão = 20)\
    -l : ComCIGAR ou SemCIGAR (padrão = SemCIGAR)
    -q : Qualidade das bases (padrão = 20)\
    -p : Porcentagem de bases cobertas (padrão = 70)\
    -c : Cobertura de reads (padrão = 20)\
    -a : Tipo de amostra analisada SP (suposto pai), M (mãe) ou P (plasma)\
    -e : Limite para erros de sequenciamento (Padrão = 10)\
    -s : Limite para considerar HOMOZIGOTO (Padrão = 80)\
    -d : Limite para considerar HETEROZIGOTO quando existem 3 ou mais (Padrão = 35)\
\n";
}

my $nome_genotipo = "Genotipos_$amostra.M$map.$cigar.Q$qual.P$por.C$cob.E$erros.S$superior.D$duvida.txt";
my $nome_report = "MCoverage_$amostra.M$map.$cigar.Q$qual.P$por.C$cob.E$erros.S$superior.D$duvida.txt";

#Usamos esse match para armazenar o início do nome de cada arquivo
$BAM =~ m/(IonXpress[\._]{1}([0-9]{3}))[\._]{1}R[\._]{1}([0-9]{4}_[0-9]{2}_[0-9]{2})((.)*).bam/;
my $nome = $1;
my $id = $2;
my $data = $3;

##Roda todos os scripts para análise

system ("/home/jaque/Desktop/Scripts_Martin/1_filtra_alinhamento.pl -b $BAM");

system ("/home/jaque/Desktop/Scripts_Martin/2_filtra_mapeamento.pl -b $BAM -m $map");

if ($cigar eq "SemCIGAR") {
    system ("/home/jaque/Desktop/Scripts_Martin/3_filtra_cobertura_microhaplotipo_sem_cigar.pl -b $BAM -m $map");

}

elsif ($cigar eq "ComCIGAR") {
    system ("/home/jaque/Desktop/Scripts_Martin/3_filtra_cobertura_microhaplotipo_com_cigar.pl -b $BAM -m $map");

}
    
system ("/home/jaque/Desktop/Scripts_Martin/4_filtra_qualidade_bases.pl -b $BAM -m $map -l $cigar -q $qual -p $por");

   
system ("/home/jaque/Desktop/Scripts_Martin/5_haplotipo_freq_micro.pl -b $BAM -m $map -l $cigar -q $qual -p $por");

system ("mkdir $nome");

system ("mv $nome.* $nome");

if (($amostra eq "SP") || ($amostra eq "M")) {
    system ("/home/jaque/Desktop/Scripts_Martin/6_extrai_genotipo.pl -b $BAM -m $map -l $cigar -q $qual -p $por -c $cob -e $erros -s $superior -d $duvida > $nome_genotipo");

}

system("/home/jaque/Desktop/Scripts_Martin/7_faz_report.pl -b $BAM -m $map -l $cigar -q $qual -p $por -c $cob -a $amostra -e $erros -s $superior -d $duvida > $nome_report");



