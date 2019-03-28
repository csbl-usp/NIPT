#!/usr/bin/perl -w

#Author: Jaqueline Wang
#MsC in Bioinformatics Graduate Program - USP

#SCRIPT TO CALCULATE THE MICROHAPLOTYPE FREQUENCIES FROM 1000 GENOMES 

#The script requires 3 parameters
#1 - Vcf file
#2 - List of SNPs that compose the microhaplotype
#3 - Population or super population initial or ALL

####################################################################################

use strict;
use Getopt::Long;

my $help = 0;
my ($vcf, $micro);
my $pop = "ALL"; 
my $super = "ALL";

GetOptions("help|h" => \$help,
	   "v=s" => \$vcf,
	   "m=s" => \$micro,
	   "p=s" => \$pop,
	   "s=s" => \$super
) or die "Failed to take the options! \n";

if ($help || !($vcf && $micro)) {die "\
This script requires two parameters. \
The output is the populational frequency of haplotypes in 1000 Genomes. \ 
To separate from population or super-population, use the parameter. \
\
Parameters: \
	-h	Show the options \
	-v	Vcf file from 1000 Genomes \
	-m	List of SNPs that compose the microhaplotype \
	-p	Population initial or ALL (do not use with -s) \
	-s	Super-population initial or ALL (do not use with -p) \
\n";
}

####################################################################################

#Armazenamos os indivíduos que foram selecionados (tanto na superpopulação quanto na população
my @sample;

if (($pop eq "ALL") && ($super eq "ALL") ){
    open (PAINEL, "Files/integrated_call_samples_v3.20130502.ALL.panel") or die "Failed to open PANEL file! \n";
    <PAINEL>;

    while (my $linha = <PAINEL>){
	my @array = split(/\t/, $linha);
	push @sample, $array[0];
    }#while
    close (PAINEL);
}#if

elsif (($pop =~ m/[A-Z]{3}/) && ($super eq "ALL")){
    open (PAINEL, "Files/integrated_call_samples_v3.20130502.ALL.panel") or die "Failed to open PANEL file! \n";
    <PAINEL>;
    while (my $linha = <PAINEL>){
	my @array = split(/\t/, $linha);
	
	if ($pop eq $array[1]){
	    push @sample, $array[0]
	}#if

    }#while
    close (PAINEL)
}#if

elsif (($super =~ m/[A-Z]{3}/) && ($pop eq "ALL")){
    open (PAINEL, "Files/integrated_call_samples_v3.20130502.ALL.panel") or die "Failed to open PANEL file! \n";
    <PAINEL>;
    while (my $linha = <PAINEL>){
	my @array = split(/\t/, $linha);
	
	if ($super eq $array[2]){
	    push @sample, $array[0]
	}#if

    }#while
    close (PAINEL)
}#if


#Abrimos o arquivo com a lista de SNPs do micro-haplótipo
open (SNPS, "$micro") or die "Failed to open MICRO file! \n";
my @snps;
while (my $entrada = <SNPS>){
    chomp ($entrada);
    push @snps, $entrada;
}
close (SNPS);

#Abrimos o arquivo VCF
$vcf =~ m/(\d*)\.(\d*)-(\d*)\.(.)*\.vcf/;
my $chrom = $1;
my $inicio = $2;
my $final = $3;

print "Chr$chrom:$inicio-$final\n";
my @micro_vcf;
my @dados;
my @ref;
my @alt;
my @id_snp;

open (TEXT, "$vcf") or die "Failed to open $vcf! \n";
my @armaz;
while (my $line = <TEXT>){
    chomp ($line);

    if ($line =~ m/\#CHROM/){
	my @amostras = split (/\t/, $line);

	foreach my $um(@amostras){
	    if ($um =~ m/[HGNA]{2}[0-9]{5}/){
		push @armaz, $um;
	    }#if
	}#foreach
    }#if

    elsif ($line =~ m/##/){
	my $string = 0;
    }#elsif

    else{
	@dados = split(/\t/, $line);
	foreach my $rs(@snps){
	    
	    if ($rs eq $dados[2]){
		push @id_snp, $dados[2];
		push @ref, $dados[3];
		push @alt, $dados[4];
		my @matriz;
		
		foreach my $GT(@dados){
		    if ($GT =~ m/\d\|\d/){
			push @matriz, $GT;
		    }#if
		}#foreach
		
		my $juntar = join ("\t", @matriz);
		push @micro_vcf, $juntar;
	    }#if
	}#foreach
    }#else
}#while	
close (TEXT);

#Calculamos as frequencias dos microhaplótipos dos indivíduos de acordo com as populações e super-populações selecionadas
my %haplotipo;
my $total = 0;

foreach my $comp(@sample){
    my $count = 0;
    while ($count < 2504){
	my @array1;
	my @array2;
	
	if ($armaz[$count] eq $comp){
	    $total += 2;
	    my $contador = 0;
	    
	    foreach my $line(@micro_vcf){
		my @dado = split(/\t/, $line);
		$dado[$count] =~ m/(\d)\|(\d)/;
		my $esquerda = $1;
		my $direita = $2;
		
		if ($alt[$contador] =~ m/[ATCG]{1},[ATCG]{1}/) {

		    $alt[$contador] =~ m/([ATCG]{1}),([ATCG]{1})/;
		    my $alt1 = $1;
		    my $alt2 = $2;
		    if ($esquerda == 0) {
			push @array1, $ref[$contador];
		    }

		    elsif ($esquerda == 1) {
			push @array1, $alt1;
		    }

		    elsif ($esquerda == 2) {
			push @array1, $alt2;
		    }

		    if ($direita == 0) {
			push @array2, $ref[$contador];
		    }

		    elsif ($direita == 1) {
			push @array2, $alt1;
		    }

		    elsif ($direita == 2) {
			push @array2, $alt2;
		    }
		}

		else {
		    if ($esquerda == 0) {
			push @array1, $ref[$contador];
		    }

		    elsif ($esquerda == 1) {
			push @array1, $alt[$contador];
		    }

		    if ($direita == 0) {
			push @array2, $ref[$contador];
		    }

		    elsif ($direita == 1) {
			push @array2, $alt[$contador];
		    }
		}

		$contador += 1;
	    }#foreach
	}#if

	else {
	    push @array1, "Nada";
	    push @array2, "Nada";
	}#else

	if (($array1[0] ne "Nada") && ($array2[0] ne "Nada")){
	    my $haplo1 = join ("", @array1);
	    my $haplo2 = join ("", @array2);
	
	    if (!exists($haplotipo{$haplo1})){
		$haplotipo{$haplo1} = 1;
	    }#if
	
	    elsif (exists($haplotipo{$haplo1})){
		$haplotipo{$haplo1} += 1;
	    }#elsif
	
	    if (!exists($haplotipo{$haplo2})){
		$haplotipo{$haplo2} = 1;
	    }#if
	    
	    elsif (exists($haplotipo{$haplo2})){
		$haplotipo{$haplo2} += 1;
	    }#elsif	
	}#if

	$count = $count + 1;
    }#while
}#foreach

#Imprimimos as frequências calculadas para cada haplótipo
foreach my $key(keys(%haplotipo)){
    my $freq = $haplotipo{$key};
    print "$key\t$freq\n";
}#foreach
