#!/usr/bin/perl -w

#Autora: Jaqueline Wang
#Mestre do Programa Interunidades em Bioinformática - USP

#SCRIPT PARA EXTRAIR OS GENÓTIPOS DO SUPOSTO PAI E DA MÃE PARA CADA MICRO-HAPLÓTIPO

#Devem ser passados 6 parâmetros:
#1 - Arquivo BAM
#2 - Qualidade do mapeamento
#3 - ComCIGAR ou SemCIGAR
#4 - Qualidade das BASES
#5 - Porcentagem das BASES cobertas
#6 - Cobertura das regiões

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
my $BAM;
my $map = 20;
my $cigar = "SemCIGAR";
my $qual = 20;
my $por = 70;
my $cob = 20;
my $erros = 10;
my $superior = 80;
my $duvida = 35;

GetOptions("help|h" => \$help,
	   "b=s" => \$BAM,
	   "m=s" => \$map,
	   "l=s" => \$cigar,
	   "q=s" => \$qual,
	   "p=s" => \$por,
	   "c=s" => \$cob,
	   "e=s" => \$erros,
	   "s=s" => \$superior,
	   "d=s" => \$duvida
    ) or die "Erro ao pegar as opções! \n";

if ($help || !($BAM)) {die "\
$0: \
Esse script recebe seis entradas, o arquivo BAM, a qualidade de mapeamento, com ou sem CIGAR, qualidade das bases, porcentagem de bases cobertas e a cobertura. A saída são os genótipos encontrados para cada micro-haplótipo a ser analisado nessa amostra.\
\
Parâmetros:\
    -h ou --help : Mostra as opções\
    -b : Arquivo BAM\
    -m : Qualidade do mapeamento dos reads (padrão = 20)\
    -l : ComCIGAR ou SemCIGAR (padrão = SemCIGAR)
    -q : Qualidade das bases (padrão = 20)\
    -p : Porcentagem de bases cobertas (padrão = 70)\
    -c : Cobertura de reads (padrão = 20)\
\
Caso queira mudar os desbalanços, pode fazer as alterações utilizando 3 variáveis:
    -e : Limite para erros de sequenciamento (Padrão = 10)\
    -s : Limite para considerar HOMOZIGOTO (Padrão = 80)\
    -d : Limite para considerar HETEROZIGOTO quando existem 3 ou mais (Padrão = 35)\
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

my $inferior = 100 - $superior;

while (my $pos = <POS>) {

    chomp ($pos);

    #Armazenamos o cromossomo e a posição inicial e final do microhaplótipo
    $pos =~ m/(chr[0-9]+):(.*)-(.*)/;
    my $chr = $1;
    my $ini = $2;
    my $fim = $3;
    
    #Abrimos o arquivo de entrada com as porcentagens relativas dos haplótipos
    open (INFILE, "$nome/$nome.haplotipos.MAIS_$por.Q$qual.M$map.$cigar.tsv") or die "Não foi possível abrir o arquivo com os haplótipos! \n";
    
    my %cobertura;
    my %porcentagem;
    
    while (my $line1 = <INFILE>) {

	chomp ($line1);
	my @dados = split(/\t/, $line1);

	if ($pos eq $dados[0]) {
	    $porcentagem{$dados[1]} = $dados[3];
	    $cobertura{$dados[1]} = $dados[2];
	}

    }
    
    close (INFILE);

    my $count10 = 0;
    my $count80 = 0;
    
    foreach my $key1(keys(%porcentagem)) {

	if ($porcentagem{$key1} > $erros) {
	    $count10 += 1;
	}

	if ($porcentagem{$key1} > $superior) {
	    $count80 += 1;
	}
    }

    print "$pos";
    my @haplos;
    
    if ($count80 == 1) {
	
	foreach my $key2(keys(%porcentagem)) {
	    
	    if (($porcentagem{$key2} > $superior) && ($cobertura{$key2} >= $cob)) {
		
		my $num = sprintf "%.2f", $porcentagem{$key2};
		my $string = "$key2\t$num%\t$cobertura{$key2}X\n";
		push @haplos, $string;
	    }

	}
	
	if (scalar @haplos == 1) {
	    print "\t";
	    print "Homozigoto\n";
	    print @haplos;
	}

	else {
	    print "\n";
	    print "Qualidade baixa! \n";


	}


    }

    elsif ($count80 == 0) {

	if ($count10 == 2) {
	    foreach my $key2(keys(%porcentagem)) {

		if (($porcentagem{$key2} >= $inferior) && ($cobertura{$key2} >= $cob)) {

		    my $num = sprintf "%.2f", $porcentagem{$key2};
		    my $string = "$key2\t$num%\t$cobertura{$key2}X\n";
		    push @haplos, $string;
		}
	    }

	    if (scalar @haplos == 2) {
		print "\t";
		print "Heterozigoto\n";
		print @haplos;
	    }

	    else {
		print "\n";
		print "Qualidade baixa! \n";
	    }
	}

	elsif ($count10 >= 3) {
	    foreach my $key2(keys(%porcentagem)) {

		if (($porcentagem{$key2} >= $duvida) && ($cobertura{$key2} >= $cob)) {

		    my $num = sprintf "%.2f", $porcentagem{$key2};
		    my $string = "$key2\t$num%\t$cobertura{$key2}X\n";
		    push @haplos, $string;
		}
	    }

	    if (scalar @haplos == 2) {
		print "\t";
		print "Heterozigoto\n";
		print @haplos;
	    }

	    else {
		print "\n";
		print "Qualidade baixa! \n";
	    }
	}

	elsif ($count10 == 1) {
	    print "\n";
	    print "Qualidade baixa! \n";
	}

	elsif ($count10 == 0) {
	    print "\n";
	    print "Qualidade baixa! Não detectou haplotipos \n";
	}
    }


    print "\n";
    
} #while (my $pos = <POS>)



close (POS);

