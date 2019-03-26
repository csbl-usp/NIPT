#!/usr/bin/perl -w

#Autora: Jaqueline Wang
#Mestre do Programa Interunidades em Bioinformática - USP

#SCRIPT PARA EXTRAIR A QUANTIDADE DE READS QUE PASSARAM NOS CONTROLES DE QUALIDADE

#Devem ser passados 7 parâmetros:
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
Esse script recebe os parâmetros utilizados para as análises e tem como saída um REPORT com as quantidades de reads que passaram em cada uma das etapas de qualidade. \
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

####################################################################################

#Usamos esse match para armazenar o início do nome de cada arquivo
$BAM =~ m/(IonXpress[\._]{1}([0-9]{3}))[\._]{1}R[\._]{1}([0-9]{4}_[0-9]{2}_[0-9]{2})((.)*).bam/;
my $nome = $1;
my $id = $2;
my $data = $3;

#Armazenamos os cromossomos e intervalos escolhidos como micro-haplótipos
open (POS, "/home/jaque/Desktop/Scripts_Martin/Arquivos/microhaplotipos.txt") or die "Não foi possível obter os micro-haplótipos! \n";

my $num1 = 1;
while (my $pos = <POS>) {

    my $micro;
    if ($num1 < 10){
	$micro = "M0".$num1;
    }

    else {
	$micro = "M".$num1;
    }
    
    chomp ($pos);

    #Armazenamos o cromossomo e a posição inicial e final do microhaplótipo    
    $pos =~ m/(chr[0-9]+):(.*)-(.*)/;
    my $chr = $1;
    my $ini = $2;
    my $fim = $3;

    system("wc $nome/$nome.MAU_ALN/$nome.$chr.$ini-$fim.MAU_ALN.tsv > report_$micro");
    
    system("wc $nome/$nome.BOM_ALN/$nome.$chr.$ini-$fim.BOM_ALN.tsv >> report_$micro");
    
    system("wc $nome/$nome.MAU_MAP.M$map/$nome.$chr.$ini-$fim.MAU_MAP.M$map.tsv >> report_$micro");
    
    system("wc $nome/$nome.BOM_MAP.M$map/$nome.$chr.$ini-$fim.BOM_MAP.M$map.tsv >> report_$micro"); 

    if ($cigar eq "SemCIGAR") {
	system("wc $nome/$nome.CIGAR_RUIM.M$map.$cigar/$nome.$chr.$ini-$fim.CIGAR_RUIM.M$map.$cigar.tsv >> report_$micro");
    }
    
    system("wc $nome/$nome.PARTE_SNPS.M$map.$cigar/$nome.$chr.$ini-$fim.PARTE_SNPS.M$map.$cigar.tsv >> report_$micro");    

    system("wc $nome/$nome.TODOS_SNPS.M$map.$cigar/$nome.$chr.$ini-$fim.TODOS_SNPS.M$map.$cigar.tsv >> report_$micro");

    system("wc $nome/$nome.MENOS_$por.Q$qual.M$map.$cigar/$nome.$chr.$ini-$fim.MENOS_$por.Q$qual.M$map.$cigar.tsv >> report_$micro");

    system("wc $nome/$nome.MAIS_$por.Q$qual.M$map.$cigar/$nome.$chr.$ini-$fim.MAIS_$por.Q$qual.M$map.$cigar.tsv >> report_$micro");

    $num1 += 1;

	
}
close (POS);


my @PAIRING;
open (POS, "/home/jaque/Desktop/Scripts_Martin/Arquivos/microhaplotipos.txt") or die "Não foi possível obter os micro-haplótipos! \n";

while (my $pos = <POS>) {

    chomp ($pos);

    open (HAPLO, "$nome/$nome.haplotipos.MAIS_$por.Q$qual.M$map.$cigar.tsv") or die "Não abriu o arquivo de haplotipos! \n";

    while (my $line1 = <HAPLO>) {
	chomp ($line1);
	if ($line1 =~ m/$pos\tDescart/g) {
	    my @dados = split(/\t/, $line1);
	    push @PAIRING, $dados[2];
	}

    } 
    close (HAPLO);
}
close (POS);


my $num2 = 1;
my @head;
my @MAU_ALN;
my @BOM_ALN;
my @MAU_MAP;
my @BOM_MAP;
my @CIGAR_RUIM;
my @PARTE_SNPS;
my @TODOS_SNPS;
my @MENOS_70;
my @MAIS_70;

while ($num2 < $num1) {

    my $micro;
    if ($num2 < 10) {
	$micro = "M0".$num2;
    }
    else {
	$micro = "M".$num2;
    }

    push @head, $micro;
    open (IN, "report_$micro") or die "Não foi possível abrir o arquivo report_$micro! \n";

    while (my $line1 = <IN>) {
	
	chomp ($line1);

	$line1 =~ m/(\s*)([0-9]+)(\s+)([0-9]+)(\s+)([0-9]+)(\s+)(.*)(\n*)/;

	my $reads = $2;
	my $arq = $8;

	$arq =~ m/$nome\/$nome\.(.*)\/$nome\.(.*)tsv/;

	my $origem = $1;

	if ($origem eq "MAU_ALN") {
	    push @MAU_ALN, $reads;
	}

	elsif ($origem eq "BOM_ALN") {
	    push @BOM_ALN, $reads;
	}

	elsif ($origem eq "MAU_MAP.M$map") {
	    push @MAU_MAP, $reads;
	}

	elsif ($origem eq "BOM_MAP.M$map") {
	    push @BOM_MAP, $reads;
	}

	elsif ($origem eq "CIGAR_RUIM.M$map.$cigar") {
	    push @CIGAR_RUIM, $reads;
	}

	elsif ($origem eq "PARTE_SNPS.M$map.$cigar") {
	    push @PARTE_SNPS, $reads;
	}

	elsif ($origem eq "TODOS_SNPS.M$map.$cigar") {
	    push @TODOS_SNPS, $reads;
	}

	elsif ($origem eq "MENOS_$por.Q$qual.M$map.$cigar") {
	    push @MENOS_70, $reads;
	}

	elsif ($origem eq "MAIS_$por.Q$qual.M$map.$cigar") {
	    push @MAIS_70, $reads;
	}
    }

    close (IN); 

    $num2 += 1;
}

my $header = join("\t", "Parameter", @head);
print "$header\n";

my $num3 = 0;
my @TOTAL;
while ($num3 < scalar @MAU_ALN) {
    my $tot = $MAU_ALN[$num3] + $BOM_ALN[$num3];
    push @TOTAL, $tot;
    $num3 += 1;
}

my $string0 = join("\t", "NUMBER_OF_READS", @TOTAL);
print "$string0\n";

my $string2 = join("\t", "AFTER_ALIGNMENT", @BOM_ALN);
print "$string2\n";

my $string4 = join("\t", "AFTER_MAPPING", @BOM_MAP);
print "$string4\n";

my @CIGAR_BOM;
my $num4 = 0;
while ($num4 < scalar @PARTE_SNPS) {
    my $tot = $PARTE_SNPS[$num4] + $TODOS_SNPS[$num4];
    push @CIGAR_BOM, $tot;
    $num4 += 1;
}

my $string = join("\t", "AFTER_CIGAR", @CIGAR_BOM);
print "$string\n";

my $string9 = join("\t", "AFTER_$por%_COB", @MAIS_70);
print "$string9\n";

my @ONE_MATCH;
my $num5 = 0;
while ($num5 < scalar @PAIRING) {
    my $tot = $MAIS_70[$num5] - $PAIRING[$num5];
    push @ONE_MATCH, $tot;
    $num5 += 1;
}

my $string11 = join("\t", "AFTER_PAIRING", @ONE_MATCH);
print "$string11\n";

my @GENOTYPE;

if (($amostra eq "M") || ($amostra eq "SP")) {

    open (POS, "/home/jaque/Desktop/Scripts_Martin/Arquivos/microhaplotipos.txt") or die "Não foi possível obter os micro-haplótipos! \n";

    while (my $pos = <POS>) {

	chomp ($pos);

	open (GENO, "Genotipos_$amostra.M$map.$cigar.Q$qual.P$por.C$cob.E$erros.S$superior.D$duvida.txt") or die "Não abriu o arquivo de genotipos! \n";
	$/ = "\n\n";
	
	while (my $line1 = <GENO>) {
	    chomp ($line1);
	    if ($line1 =~ m/$pos/g) {
		my @array = split(/\n/, $line1);
		my $sum = 0;

		foreach my $key (@array) {

		    if ($key =~ m/[ATCG]+/g) {

			my ($haplo, $porce, $cober) = split(/\t/, $key);
			$cober =~ m/([0-9]+)X/;
			my $cob = $1;
			$sum += $cob;
			
		    }
		}

		push @GENOTYPE, $sum;
	    }
	} 
	close (GENO);

	$/ = "\n";
    }
    close (POS);
    my $string12 = join("\t", "AFTER_GENOTYPE", @GENOTYPE);
    print "$string12\n";
    
}
system ("rm report_M*");
