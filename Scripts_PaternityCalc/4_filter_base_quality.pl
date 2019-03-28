#!/usr/bin/perl -w

#Author: Jaqueline Wang
#MsC in Bioinformatics Graduate Program - USP

#SCRIPT TO EXTRACT THE HAPLOTYPES FROM THE READS WITH GOOD MAPPING AND COVERING ALL THE SNPS WITHIN THE MICROHAPLOTYPE 

#The script receives 7 parameters
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
The outputs are two files. \
	MORE_XX - contém os haplótipos com mais de XX% de bases cobertas. \
	LESS_XX - contém os haplótipos com mais de XX% de bases cobertas. \ 
\
Required parameters: \
	-h	Mostra as opções \
	-b	Bam file \
	-t	Trio number \
	-a	Type of sample AF (alleged father), M (mother) ou P (plasma) \ 
\
Other parameters: \
	-m	Mapping quality of reads (default = 20) \
	-l	Uses or not the CIGAR info. If NOT (NotCIGAR), consider only (mis)matches. Other option is YesCIGAR (default = NotCIGAR) \
	-q	Bases quality (default = 20) \
	-p	Percentage of covered bases (default = 70) \
\n";
}

####################################################################################

my $nome = "Trio$trio"."_Sample$amostra";

my %qualidade;
open (QUAL, "Files/phred.txt") or die "Failed to open file with bases quality library! \n";

while (my $line1 = <QUAL>) {
    chomp ($line1);
    my ($qual1, $qual2) = split (/\t/, $line1);

    if (!exists $qualidade{$qual1}) {
	$qualidade{$qual1} = $qual2;
    }
} #while (my $line1 = <QUAL>)

close (QUAL);

open (POS, "Files/Microhaplotypes.txt") or die "Failed to open microhaplotypes! \n";

while (my $pos = <POS>) {
    chomp ($pos);
    
    #Usamos esse match para armazenar o cromossomo e a posição inicial e final do micro-haplótipo
    $pos =~ m/(chr[0-9]+):(.*)-(.*)/;
    my $chr = $1;
    my $ini = $2;
    my $fim = $3;

    #Abrimos o arquivo com TODOS os SNPs cobertos
    open (INFILE, "$nome.TODOS_SNPS.M$map.$cigar/$nome.$chr.$ini-$fim.TODOS_SNPS.M$map.$cigar.tsv") or die "Failed to open TODOS_SNPS file! \n";
    
    my @micro;
    my $soma = 0;
     
    while (my $line2 = <INFILE>) {
	chomp ($line2);
	$soma = $soma + 1;
	my ($haplo1, $qual) = split(/\t/, $line2);
	
	my @bases = $haplo1 =~ m/[-ATCG]/g;
	my @dados = split(//, $qual);
	
	my $tamanho = scalar(@bases);
	my $count = 0;
	my @new_bases;
	my @new_qual1;
	
	while ($count < $tamanho) {
	    foreach my $key1(keys(%qualidade)) {
		if (($dados[$count] eq $key1) && ($dados[$count] == $key1)){
		    
		    if ($qualidade{$key1} >= $score) {
			push @new_bases, $bases[$count];
			push @new_qual1, $qualidade{$key1};
		    }
		    
		    elsif ($qualidade{$key1} < $score) {
			my $erro = "-";
			push @new_bases, $erro;
			push @new_qual1, "0";
		    }
		}
	    }
	    
	    $count = $count + 1;
	}

	my $haplo = join("",@new_bases);
	my $qual_info1 = join(",", @new_qual1);

	my $concat = join("\t", $haplo, $qual_info1);
	push @micro, $concat;
    } #while (my $line2 = <INFILE>)
    
    close (INFILE);

    #Arquivo com as posições dos SNPs
    open (BED, "Files/SNPs.bed") or die "Failed to open BED file! \n";
    
    my @SNP_id;
    
    while (my $line3 = <BED>) {
	chomp ($line3);
	my @info = split(/\t/, $line3);
	
	if (($info[0] eq $chr) && ($info[1] >= $ini) && ($info[1] <= $fim)){
	    push @SNP_id, $info[3];
	} #if (($info[0] eq $chr) && ($info[1] >= $ini) && ($info[1] <= $fim))
    } #while (my $line3 = <BED>)
    
    close (BED); 

    #Abrimos o arquivo com parte dos SNPs cobertos
    open (INFILE1, "$nome.PARTE_SNPS.M$map.$cigar/$nome.$chr.$ini-$fim.PARTE_SNPS.M$map.$cigar.tsv") or die "Failed to open PARTE_SNPS file! \n";   
    
    while (my $line4 = <INFILE1>) {
	$soma = $soma + 1;
	chomp ($line4);
	my ($haplo2, $qual2, $SNP) = split(/\t/, $line4);
	
	my @id = split (/ /, $SNP);
	my @bases = $haplo2=~ m/[ATCG-]/g;
	my @dados = split(//, $qual2);
	
	my $count1 = 0;
	my $count2 = 0;
	my @new_bases;
	my @new_qual1;
	
	my $tam = scalar @SNP_id;
	
	while ($count1 <= $tam) {
	    
	    if ($id[$count2] eq $SNP_id[$count1]) {
		foreach my $key1(keys(%qualidade)) {
		    
		    if (($dados[$count2] eq $key1) && ($dados[$count2] == $key1)) {
			
			if ($qualidade{$key1} >= $score) {
			    push @new_bases, $bases[$count2];
			    push @new_qual1, $qualidade{$key1};
			}
			
			elsif ($qualidade{$key1} < $score) {
			    my $erro = "-";
			    push @new_bases, $erro;
			    push @new_qual1, "0";
			}
		    }
		}
		$count1 = $count1 + 1;
		$count2 = $count2 + 1;

	    } # if
	    
	    elsif ($id[$count2] ne $SNP_id[$count1]) {
		push @new_bases, "-";
		push @new_qual1, "0";
		$count1 = $count1 + 1;
	    } # elsif
	} # while ($count1 <= $tam)
	
	my $haplo1 = join("", @new_bases);
	my $qual1 = join(",", @new_qual1);
	my $concat = join("\t", $haplo1, $qual1);
	
	push @micro, $concat;
    }
    close (INFILE1);

    #Arquivo com MORE
    open (OUTMAIS, ">$nome.$chr.$ini-$fim.MORE_$por.Q$score.M$map.$cigar.tsv") or die "Failed to open MORE file! \n";

    #Arquivo com LESS
    open (OUTMENOS, ">$nome.$chr.$ini-$fim.LESS_$por.Q$score.M$map.$cigar.tsv") or die "Failed to open LESS file! \n";    

    my @mais_70;
    my @mais_100;
    my $descartados = 0;
    
    foreach my $key2(@micro) {
	my ($haplo, $qual1) = split(/\t/, $key2);
	my @matches = $haplo =~ m/[-]/g;
	my $indef = scalar(@matches);
	my $ref = scalar(@SNP_id);

	my $aux = (100 - $por) * 0.01;

	if ($indef < ($ref * $aux)) {
	    push @mais_70, $haplo;
	    print OUTMAIS "$key2\n";

	    if ($indef == 0) {
		push @mais_100, $haplo;
	    }
	}
	    
	else {
	    print OUTMENOS "$key2\n";
	    $descartados += 1;
	}
    }
    close (OUTMAIS);
    close (OUTMENOS);
} #while (my $pos = <POS>)

close (POS);

system ("mkdir $nome.MORE_$por.Q$score.M$map.$cigar $nome.LESS_$por.Q$score.M$map.$cigar");
system ("mv $nome.*.MORE_$por.Q$score.M$map.$cigar.tsv $nome.MORE_$por.Q$score.M$map.$cigar");
system ("mv $nome.*.LESS_$por.Q$score.M$map.$cigar.tsv $nome.LESS_$por.Q$score.M$map.$cigar");
