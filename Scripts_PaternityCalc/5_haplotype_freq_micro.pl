#!/usr/bin/perl -w

#Author: Jaqueline Wang
#MsC in Bioinformatics Graduate Program - USP

#SCRIPT TO SEPARATE READS THAT HAD PASSED THE ALIGNMENT, MAPPING, BASE QUALITY, SNPS COVERAGE AND PERCENTAGE OF COVERED SNPS.

#the script receives 5 parameters
#1 - Bam file
#2 - Mapping quality
#3 - YesCIGAR or NotCIGAR
#4 - Bases quality
#5 - Percentage of covered bases

####################################################################################

use strict;
use Getopt::Long;

my $help = 0;
my $BAM;
my $map = 20;
my $cigar = "NotCIGAR";
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
This script receives five inputs, the bam file, mapping quality, YesCIGAR or NotCIGAR, bases quality and percentage of covered bases. \
The output is a file with the haplotypes found in each microhaplotype.\
\
Parameters:\
    -h ou --help : Mostra as opções\
    -b : Bam file\
    -m : Mapping quality of reads (default = 20)\
    -l : Uses or not the CIGAR info. If NOT (NotCIGAR), consider only (mis)matches. Other option is YesCIGAR (default = NotCIGAR) \
    -q : Bases quality (default = 20)\
    -p : Percentage of covered bases (default = 70)\
\n";
}

####################################################################################

#Usamos esse match para armazenar o início do nome de cada arquivo
$BAM =~ m/(IonXpress[\._]{1}([0-9]{3}))[\._]{1}R[\._]{1}([0-9]{4}_[0-9]{2}_[0-9]{2})((.)*).bam/;
my $nome = $1;
my $id = $2;
my $data = $3;

#Armazenamos os cromossomos e intervalos escolhidos como micro-haplótipos
open (POS, "Files/microhaplotipos.txt") or die "Failed to open file with microhaplotypes! \n";

#Abrimos o arquivo de saída dos haplótipos encontrados em cada micro-haplótipo
open (OUT, ">$nome.haplotipos.MORE_$por.Q$score.M$map.$cigar.tsv") or die "Failed to open HAPLOTYPES file! \n";

while (my $pos = <POS>) {

    chomp ($pos);
    
    #Armazenamos o cromossomo e a posição inicial e final do microhaplótipo
    $pos =~ m/(chr[0-9]+):(.*)-(.*)/;
    my $chr = $1;
    my $ini = $2;
    my $fim = $3;

    
    open (INFILE, "$nome.MORE_$por.Q$score.M$map.$cigar/$nome.$chr.$ini-$fim.MORE_$por.Q$score.M$map.$cigar.tsv") or die "Failed to open MORE file! \n";

 
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

