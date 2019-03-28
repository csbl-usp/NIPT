#!/usr/bin/perl -w

#Author: Jaqueline Wang
#MsC in Bioinformatics Graduate Program - USP

#SCRIPT TO SEPARATE READS THAT HAD PASSED THE ALIGNMENT, MAPPING, BASE QUALITY, SNPS COVERAGE AND PERCENTAGE OF COVERED SNPS.

#the script receives 7 parameters
#1 - Bam file
#2 - Mapping quality
#3 - YesCIGAR or NotCIGAR
#4 - Bases quality
#5 - Percentage of covered bases
#6 - Trio number
#7 - Sample type

####################################################################################

use strict;
use Getopt::Long;

my $help = 0;
my $BAM;
my $map = 20;
my $cigar = "NotCIGAR";
my $score = 20;
my $por = 70;
my $trio;
my $amostra;

GetOptions("help|h" => \$help,
	   "b=s" => \$BAM,
	   "m=s" => \$map,
	   "l=s" => \$cigar,
	   "q=s" => \$score,
	   "p=s" => \$por,
	   "t=s" => \$trio,
	   "a=s" => \$amostra
    ) or die "Failed to take the options! \n";

if ($help || !($BAM && $trio && $amostra)) {die "\
This script requires three inputs, the bam file, the trio number and the sample type. \
Other parameters have default values, but can be changed. \
The output is a file with the haplotypes found in each microhaplotype.\
\
Required parameters:\
	-b	Bam file\
	-t	Trio number \
	-a	Type of sample AF (alleged father), M (mother) ou P (plasma) \ 
\
Other parameters: \
	-h	Show the options \
	-m	Mapping quality of reads (default = 20) \
	-l	Uses or not the CIGAR info. If NOT (NotCIGAR), consider only (mis)matches. Other option is YesCIGAR (default = NotCIGAR) \
	-q	Bases quality (default = 20) \
	-p	Percentage of covered bases (default = 70) \
\n";
}

####################################################################################

my $nome = "Trio$trio"."_Sample_$amostra";

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
    print OUT "$chr:$ini-$fim\tDiscard\t$descartados\n";
    print OUT "\n";
    
} #while (my $pos = <POS>)


close (OUT);
close (POS);

