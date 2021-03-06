#!/usr/bin/perl -w

# Author: Jaqueline Wang
# MsC in Bioinformatics Graduate Program - USP

# SCRIPT TO SEPARATE READS THAT COVER ALL THE SNPS CONSIDERING THE CIGAR

# The script receives 3 parameters
# 1 - Mapping quality
# 2 - Trio number
# 3 - Sample type

####################################################################################

use strict;
use Getopt::Long;

my $help = 0;
my $map = 20;
my $trio;
my $amostra;

GetOptions("help|h" => \$help,
	   "m=s" => \$map,
	   "t=s" => \$trio,
	   "a=s" => \$amostra
    ) or die "Failed to take the options! \n";

if ($help || !($trio && $amostra)) {die "\
This script requires two inputs, the trio number and the sample type. \
Other parameters have default values, but can be changed. \
The outputs are two files. \
	TODOS_SNPS - contain reads with ALL the SNPs covered. \
	PARTE_SNPS - contain reads with PART of the SNPs covered. \
\
Required parameters: \
	-t	Trio number \
	-a	Type of sample AF (alleged father), M (mother) ou P (plasma) \
\
Other parameters: \
	-h	Show the options \
	-m	Mapping quality of reads (default = 20) \
\n";
}

####################################################################################

my $nome = "Trio$trio"."_Sample$amostra";

#Armazenamos os cromossomos e intervalos escolhidos como micro-haplótipos
open (POS, "Files/Microhaplotypes.txt") or die "Failed to obtain the microhaplotypes! \n";

while (my $pos = <POS>) {
    chomp ($pos);

    #Armazenamos o cromossomo e a posição inicial e final do microhaplótipo
    $pos =~ m/(chr[0-9]+):(.*)-(.*)/;
    my $chr = $1;
    my $ini = $2;
    my $fim = $3;

    #Abrimos o arquivo com as posições dos SNPs
    open (BED, "Files/SNPs.bed") or die "Failed to open BED file! \n";

    my @SNP_pos;
    my @SNP_id;
    my $n_snp = 0;

    while (my $line1 = <BED>) {
	chomp ($line1);
	my @info = split(/\t/, $line1);

	if (($info[0] eq $chr) && ($info[1] >= $ini) && ($info[1] <= $fim)){
	    push @SNP_pos, $info[1];
	    push @SNP_id, $info[3];
	    $n_snp += 1;
	} #if (($info[0] eq $chr) && ($info[1] >= $ini) && ($info[1] <= $fim))


    } #while (my $line1 = <BED>)

    close (BED);

    #Abrimos o arquivo com BOM mapeamento
    open (INFILE, "$nome.BOM_MAP.M$map/$nome.$chr.$ini-$fim.BOM_MAP.M$map.tsv") or die "Failed to open BOM_MAP file! \n";

    my @haplo_bom;
    my @haplo_ruim;

    while (my $line2 = <INFILE>) {

	my @bam_info = split (/\t/, $line2);
	my $qname = $bam_info[0];
	my $flag = $bam_info[1];
	my $rname = $bam_info[2];
	my $pos0 = $bam_info[3];
	my $mapq = $bam_info[4];
	my $cigar = $bam_info[5];
	my $seq = $bam_info[9];
	my $qual = $bam_info[10];

       	my $pos1 = 0;
	my $seq1;
	my $qual1;

	my @CIGAR = $cigar =~ m/[0-9]*[MIDNSHP]{1}/g;

	foreach my $c_string(@CIGAR){
	    $c_string =~ m/([0-9]*)([MIDNSHP]{1})/;
	    my $number = $1;
	    my $MIDNSHP = $2;

	    if ($MIDNSHP eq "M"){
		my $seq2 = substr ($seq,$pos1,$number);
		my $qual2 = substr ($qual,$pos1,$number);
		$seq1 = $seq1.$seq2;
		$qual1 = $qual1.$qual2;
		$pos1 = $pos1 + $number;
	    }
		
	    elsif ($MIDNSHP eq "I"){
		my $seq2 = "";
		my $qual2 = "";
		$seq1 = $seq1.$seq2;
		$qual1 = $qual1.$qual2;
		$pos1 = $pos1 + $number;
	    }
		
	    elsif ($MIDNSHP eq "D"){	    
		my $seq2 = "-" x $number;
		my $qual2 = "!" x $number;
		$seq1 = $seq1.$seq2;
		$qual1 = $qual1.$qual2;
	    }
	    
	    elsif ($MIDNSHP eq "N"){
		my $seq2 = "-" x $number;
		my $qual2 = "!" x $number;
		$seq1 = $seq1.$seq2;
		$qual1 = $qual1.$qual2;
		$pos1 = $pos1;
	    }
	    
	    elsif ($MIDNSHP eq "S"){	    
		my $seq2 = substr($seq, $pos1, $number);
		my $qual2 = substr($qual, $pos1, $number);
		$seq1 = $seq1;
		$qual1 = $qual1;
		$pos1 = $pos1 + $number;
	    }
	    
	    elsif ($MIDNSHP eq "H"){
		my $seq2 = "";
		my $qual2 = "";
		$seq1 = $seq1.$seq2;
		$qual1 = $qual1.$qual2;
		$pos1 = $pos1;
	    }
	    
	    elsif ($MIDNSHP eq "P"){
		my $seq2 = "";
		my $qual2 = "";
		$seq1 = $seq1.$seq2;
		$qual1 = $qual1.$qual2;
		$pos1 = $pos1;
	    }
	} #foreach my $c_string(@CIGAR)

	my @hapl_seq;
	my @hapl_qual;
	my @hapl_snp;
	my $count1 = 0;
	my $count2 = 0;
	foreach my $snp_pos(@SNP_pos){
	    
	    my $base;
	    my $quab;
	    
	    if ($snp_pos >= $pos0){		    
		my $num = $snp_pos - $pos0;
		$base = substr ($seq1, $num, 1);
		$quab = substr ($qual1, $num, 1);

		if ($base =~ m/[ATCG-]/){
		    push @hapl_seq, $base;
		    push @hapl_qual, $quab;
		    push @hapl_snp, $SNP_id[$count1];
		    $count2 = $count2 + 1;
		}
		$count1 = $count1 + 1;
	    } #if ($snp_pos >= $pos0)

	    elsif ($snp_pos < $pos0) {
		$count1 = $count1 + 1;
	    } #elsif ($snp_pos < $pos0)
	    
	} #foreach my $snp_pos(@SNP_pos)

	my $haplo1 = join("", @hapl_seq);
	my $haplo2 = join("", @hapl_qual);
	my $haplo3 = join(" ", @hapl_snp);
	
	if ($count2 == $n_snp) {
	    my $teste = join ("\t", $haplo1, $haplo2);
	    push @haplo_bom, $teste;
	} #if (length($haplo1) == %n_snp

	if ($count2 != $n_snp) {
	    my $teste = join ("\t", $haplo1, $haplo2, $haplo3);
	    push @haplo_ruim, $teste;
	} #if (length($haplo1) != $n_snp)

    } #while (my $line = <INFILE>)

    close (INFILE);

    #Arquivo de reads com TODOS os SNPs cobertos
    open (OUTBOM, ">$nome.$chr.$ini-$fim.TODOS_SNPS.M$map.YesCIGAR.tsv") or die "Failed to open TODOS_SNPS file! \n";

    foreach my $g_string(@haplo_bom) {
	print OUTBOM $g_string, "\n";
    } #foreach my $g_string(@haplo_bom)

    close (OUTBOM);

    #Arquivo de reads com PARTE dos SNPs cobertos
    open (OUTRUIM, ">$nome.$chr.$ini-$fim.PARTE_SNPS.M$map.YesCIGAR.tsv") or die "Failed to open PARTE_SNPS file! \n";

    foreach my $b_string(@haplo_ruim) {
	print OUTRUIM $b_string, "\n";
    } #foreach my $b_string(@haplo_ruim)
    close (OUTRUIM);
} #while (my $pos = <POS>)

close (POS);
system ("mkdir $nome.PARTE_SNPS.M$map.YesCIGAR $nome.TODOS_SNPS.M$map.YesCIGAR");
system ("mv $nome.*.PARTE_SNPS.M$map.YesCIGAR.tsv $nome.PARTE_SNPS.M$map.YesCIGAR");
system ("mv $nome.*.TODOS_SNPS.M$map.YesCIGAR.tsv $nome.TODOS_SNPS.M$map.YesCIGAR");
