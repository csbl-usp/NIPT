#!/usr/bin/perl -w

#Autora: Jaqueline Wang
#Mestre do Programa Interunidades em Bioinformática - USP

#SCRIPT PARA SEPARAR OS READS QUE TENHAM PASSADO PELO ALINHAMENTO, MAPEAMENTO, QUALIDADE DAS BASES, COBERTURA DOS SNPS E PORCENTAGEM DE SNPS COBERTOS.

#Devem ser passados 5 parâmetros:
#1 - Arquivo BAM
#2 - Qualidade do mapeamento
#3 - ComCIGAR ou SemCIGAR
#4 - Qualidade das BASES
#5 - Porcentagem de BASES cobertas

#Endereços importantes (MARTIN)
#Micro-haplótipos: /home/jaque/Desktop/Scripts_Martin/Arquivos/microhaplotipos.txt

####################################################################################

use strict;
use Getopt::Long;

my $help = 0;
my $BAM;
my $map = 20;
my $cigar = "SemCIGAR";
my $score = 20;
my $por = 70;

GetOptions("help|h" => \$help,
	   "b=s" => \$BAM,
	   "m=s" => \$map,
	   "l=s" => \$cigar,
	   "q=s" => \$score,
	   "p=s" => \$por
    ) or die "Erro ao pegar as opções! \n";

if ($help || !($BAM)) {die "\
$0: \
Esse script recebe cinco entradas, o arquivo BAM, a qualidade de mapeamento, com ou sem CIGAR, qualidade das bases e porcentagem das bases cobertas. A saída é um arquivo com os haplótipos encontrados para cada micro-haplótipo.\
\
Parâmetros:\
    -h ou --help : Mostra as opções\
    -b : Arquivo BAM\
    -m : Qualidade do mapeamento dos reads (padrão = 20)\
    -l : ComCIGAR ou SemCIGAR (padrão = SemCIGAR)
    -q : Qualidade das bases (padrão = 20)\
    -p : Porcentagem de bases cobertas (padrão = 70)\
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

#Abrimos o arquivo de saída dos haplótipos encontrados em cada micro-haplótipo
open (OUT, ">$nome.haplotipos.MAIS_$por.Q$score.M$map.$cigar.tsv") or die "Não foi possível abrir o arquivo de saída! \n";

while (my $pos = <POS>) {

    chomp ($pos);
    
    #Armazenamos o cromossomo e a posição inicial e final do microhaplótipo
    $pos =~ m/(chr[0-9]+):(.*)-(.*)/;
    my $chr = $1;
    my $ini = $2;
    my $fim = $3;

    
    open (INFILE, "$nome.MAIS_$por.Q$score.M$map.$cigar/$nome.$chr.$ini-$fim.MAIS_$por.Q$score.M$map.$cigar.tsv") or die "Não foi possível abrir o arquivo com os haplótipos! \n";

 
    my %haplotipos;
    my @parcial;
    my $soma = 0;
    
    while (my $line1 = <INFILE>) {
	chomp ($line1);
	my ($haplo, $qual) = split (/\t/, $line1);

	my @indef = $haplo =~ m/[-]/g;
	my $cont = scalar(@indef);
	
	if ($cont == 0) {

	    if (!exists($haplotipos{$haplo})) {
		$haplotipos{$haplo} = 1;

	    }

	    elsif (exists($haplotipos{$haplo})) {
		$haplotipos{$haplo} += 1;

	    }

	}

	else {
	    push @parcial, $haplo;

	}

	$soma += 1;
	
    }
    
    close (INFILE);

    my $descartados = 0;
    
    foreach my $key1(@parcial) {
	
	$key1 =~ s/-/\./g;

	my $cont = 0;
	my $haplo;
	foreach my $key2(keys(%haplotipos)) {
	    if ($key2 =~ m/$key1/g) {
		$cont += 1;

	    }
	}

	if ($cont == 1) {
	    foreach my $key3(keys(%haplotipos)) {
		if ($key3 =~ m/$key1/g) {
		    $haplotipos{$key3} += 1;
		}
	    }
	}

	else {
	    $descartados += 1;
	}
    }
    my $div = $soma - $descartados; 
    foreach my $key(keys(%haplotipos)) {
	my $porcent = ($haplotipos{$key}/$div) * 100;
	print OUT "$chr:$ini-$fim\t$key\t$haplotipos{$key}\t$porcent\n"

    }
    print OUT "$chr:$ini-$fim\tDescart\t$descartados\n";
    print OUT "\n";
    
} #while (my $pos = <POS>)


close (OUT);
close (POS);

