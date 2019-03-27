#!/usr/bin/perl -w

#Author: Jaqueline Wang
#MsC in Bioinformatics Graduate Program - USP

#SCRIPT TO EXTRACT THE READS QUALITY THAT PASSED THE QUALITY CONTROL

#The script receives 7 parameters
#1 - Bam file
#2 - Mapping quality
#3 - YesCIGAR or NotCIGAR
#4 - Bases quality
#5 - Percentage of covered bases
#6 - Coverage of regions
#7 - Sample type

#Inbalance parameters
#1 - Superior
#2 - Errors
#3 - Doubt between three

####################################################################################

use strict;
use Getopt::Long;

my $help = 0;
my $BAM;
my $map = 20;
my $cigar = "NotCIGAR";
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
This script receives the parameters used for analysis and make a REPORT with the number of reads thas passes each quality step.\
\
Parameters:\
    -h ou --help : Show the options\
    -b : Bam file
    -m : Mapping quality of read (default = 20)\
    -l : YesCIGAR or NotCIGAR (default = NotCIGAR) \
    -q : Bases quality (default = 20)\
    -p : Percentage of covered bases (default = 70)\
    -c : Coverage (default = 20)\
    -a : Type of sample AF (alleged father), M (mother) or P (plasma)\
    -e : Limit for sequencing errors (default = 10)\
    -s : Limit to consider HOMOZYGOUS (default = 80)\
    -d : Limit to consider HETEROZYGOUS when there are more than 3 possibilities (default = 35)\
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

    system("wc $nome/$nome.LESS_$por.Q$qual.M$map.$cigar/$nome.$chr.$ini-$fim.LESS_$por.Q$qual.M$map.$cigar.tsv >> report_$micro");

    system("wc $nome/$nome.MORE_$por.Q$qual.M$map.$cigar/$nome.$chr.$ini-$fim.MORE_$por.Q$qual.M$map.$cigar.tsv >> report_$micro");

    $num1 += 1;

	
}
close (POS);


my @PAIRING;
open (POS, "Files/microhaplotypes.txt") or die "Não foi possível obter os micro-haplótipos! \n";

while (my $pos = <POS>) {

    chomp ($pos);

    open (HAPLO, "$nome/$nome.haplotipos.MORE_$por.Q$qual.M$map.$cigar.tsv") or die "Failed to open HAPLOTYPES file! \n";

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
    open (IN, "report_$micro") or die "Failed to open report_$micro file! \n";

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

	elsif ($origem eq "LESS_$por.Q$qual.M$map.$cigar") {
	    push @MENOS_70, $reads;
	}

	elsif ($origem eq "MORE_$por.Q$qual.M$map.$cigar") {
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

    open (POS, "Files/microhaplotypes.txt") or die "Failed to obtain the microhaplotypes! \n";

    while (my $pos = <POS>) {

	chomp ($pos);

	open (GENO, "Genotypes_$amostra.M$map.$cigar.Q$qual.P$por.C$cob.E$erros.S$superior.D$duvida.txt") or die "Failed to open GENOTYPES files! \n";
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
