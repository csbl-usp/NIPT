#!/usr/bin/perl -w

# Author: Jaqueline Wang
# MsC in Bioinformatics Graduate Program - USP

# SCRIPT TO EXTRACT THE GENOTYPES OF ALLEGED FATHER AND MOTHER FOR EACH MICROHAPLOTYPE

# The script receives 7 parameters
# 1 - Mapping quality
# 2 - YesCIGAR or NotCIGAR
# 3 - Bases quality
# 4 - Percentage of covered bases
# 5 - Regions coverage
# 6 - Trio number
# 7 - Sample type

# Inbalance parameters
# 1 - Superior
# 2 - Errors
# 3 - Doubt between three

####################################################################################

use strict;
use Getopt::Long;

my $help = 0;
my $map = 20;
my $cigar = "NotCIGAR";
my $qual = 20;
my $por = 70;
my $cob = 20;
my $erros = 10;
my $superior = 80;
my $duvida = 35;
my $trio;
my $amostra;

GetOptions("help|h" => \$help,
	   "m=s" => \$map,
	   "l=s" => \$cigar,
	   "q=s" => \$qual,
	   "p=s" => \$por,
	   "c=s" => \$cob,
	   "e=s" => \$erros,
	   "s=s" => \$superior,
	   "d=s" => \$duvida,
	   "t=s" => \$trio,
	   "a=s" => \$amostra
    ) or die "Failed to take the options! \n";

if ($help || !($trio && $amostra)) {die "\
This script requires two inputs, the trio number and the sample type. \
Other parameters have default values, but can be changed. \
The output is the genotype of each microhaplotype analyzed.\
\
Required parameters: \
	-t	Trio number \
	-a	Type of sample AF (alleged father) or M (mother) \
\
Other parameters: \
	-h 	Show the options \
	-m	Mapping quality of reads (default = 20) \
	-l	Uses or not the CIGAR info. If NOT (NotCIGAR), consider only (mis)matches. Other option is YesCIGAR (default = NotCIGAR) \
	-q	Bases quality (default = 20) \
	-p	Percentage of covered bases (default = 70) \
	-c	Coverage (default = 20) \
\
Inbalance parameters:
	-e	Limit for sequencing errors (default = 10)\
	-s	Limit to consider HOMOZYGOUS (default = 80)\
	-d	Limit to consider HETEROZYGOUS when ther are more than 3 possibilities (default = 35)\
\n";
}

####################################################################################

my $nome = "Trio$trio"."_Sample$amostra";

#Armazenamos os cromossomos e intervalos escolhidos como micro-haplótipos
open (POS, "Files/Microhaplotypes.txt") or die "Failed to obtain the microhaplotypes! \n";

my $inferior = 100 - $superior;

while (my $pos = <POS>) {

    chomp ($pos);

    #Armazenamos o cromossomo e a posição inicial e final do microhaplótipo
    $pos =~ m/(chr[0-9]+):(.*)-(.*)/;
    my $chr = $1;
    my $ini = $2;
    my $fim = $3;
    
    #Abrimos o arquivo de entrada com as porcentagens relativas dos haplótipos
    open (INFILE, "BAM/Trio$trio/$nome/$nome.haplotipos.MORE_$por.Q$qual.M$map.$cigar.tsv") or die "Failed to open the HAPLOTYPES file! \n";
    
    my %cobertura;
    my %porcentagem;
    
    while (my $line1 = <INFILE>) {
    
    	next unless $line1 !~ m/^\n/;
	chomp ($line1);
	my @dados = split(/\t/, $line1);

	if ($pos eq $dados[0]) {
	    if ($dados[1] ne "Discard") {
		$porcentagem{$dados[1]} = $dados[3];
		$cobertura{$dados[1]} = $dados[2];
	    }
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
	    print "Homozygous\n";
	    print @haplos;
	}

	else {
	    print "\n";
	    print "Low quality! \n";
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
		print "Heterozygous\n";
		print @haplos;
	    }

	    else {
		print "\n";
		print "Low quality! \n";
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
		print "Heterozygous\n";
		print @haplos;
	    }

	    else {
		print "\n";
		print "Low quality! \n";
	    }
	}

	elsif ($count10 == 1) {
	    print "\n";
	    print "Low quality! \n";
	}

	elsif ($count10 == 0) {
	    print "\n";
	    print "Low quality! No haplotypes. \n";
	}
    }
    print "\n";
} #while (my $pos = <POS>)

close (POS);

