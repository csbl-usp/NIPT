#!/usr/bin/perl -w

#Autora : Jaqueline Wang
#Mestre do Programa Interunidades em Bioinformática - USP

#SCRIPT PARA SEPARAR OS READS QUE COBREM TODOS OS SNPS SEM CIGAR

#DeveM ser passados 2 parâmetros:
#1 - Arquivo BAM
#2 - Qualidade do MAPEAMENTO

#Endereços importantes (MARTIN)
#Micro-haplótipos: /home/jaque/Desktop/Scripts_Martin/Arquivos/microhaplotipos.txt
#Arquivo BED: /home/jaque/Desktop/Scripts_Martin/Arquivos/SNPs.bed

####################################################################################

use strict;
use Getopt::Long;

my $help = 0;
my $BAM;
my $map = 20;

GetOptions("help|h" => \$help,
	   "b=s" => \$BAM,
	   "m=s" => \$map
    ) or die "Erro ao pegar as opções! \n";

if ($help || !($BAM)) {die "\
$0: \
Esse script recebe uma entrada, o arquivo BAM. A saída são três arquivos: \
CIGAR_RUIM - contém os reads cujo CIGAR contém inserções, deleções, etc. \
TODOS_SNPS - contém os redas que possuem TODOS os SNPs cobertos. \
PARTE_SNPS - contém os reads que possuem PARTE dos SNPs cobertos. \
\
Parâmetros:\
    -h ou --help : Mostra as opções\
    -b : Arquivo BAM\
    -m : Qualidade do mapeamento dos reads (padrão = 20)\
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

while (my $pos = <POS>) {
    
    chomp ($pos);

    #Armazenamos o cromossomo e a posição inicial e final do microhaplótipo
    $pos =~ m/(chr[0-9]+):(.*)-(.*)/;
    my $chr = $1;
    my $ini = $2;
    my $fim = $3;

    #Abrimos o arquivo com as posições dos SNPs
    open (BED, "/home/jaque/Desktop/Scripts_Martin/Arquivos/SNPs.bed") or die "Não abriu o BED! \n";

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
    open (INFILE, "$nome.BOM_MAP.M$map/$nome.$chr.$ini-$fim.BOM_MAP.M$map.tsv") or die "Não foi possível abrir o arquivo com os BOMs MAPEAMENTOS! \n";

    #Arquivo com CIGAR ruim
    open (CIGAR, ">$nome.$chr.$ini-$fim.CIGAR_RUIM.M$map.SemCIGAR.tsv") or die "Não abriu o arquivo dos reads com CIGAR ruim! \n";
    my @haplo_bom;
    my @haplo_ruim;

    my $cont_ruim = 0;
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

	if ($cigar =~ m/[SHIDNP]{1}/g) {
	    $cont_ruim += 1;
	    print CIGAR "$line2";
	}

	else {
	    
	    $cigar =~ m/([0-9]*)M/g;
	    my $number = $1;
	    my $seq2 = substr ($seq,$pos1,$number);
	    my $qual2 = substr ($qual,$pos1,$number);
	    $seq1 = $seq1.$seq2;
	    $qual1 = $qual1.$qual2;
	    $pos1 = $pos1 + $number;
	   
 

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

	}
    } #while (my $line = <INFILE>)

    close (INFILE);
    close (CIGAR);

    #Arquivo de reads com TODOS os SNPs cobertos
    open (OUTBOM, ">$nome.$chr.$ini-$fim.TODOS_SNPS.M$map.SemCIGAR.tsv") or die "Não abriu o arquivo com todos os SNPs! \n";

    foreach my $g_string(@haplo_bom) {
	print OUTBOM $g_string, "\n";
	
    } #foreach my $g_string(@haplo_bom)

    close (OUTBOM);

    #Arquivo de reads com PARTE dos SNPs cobertos
    open (OUTRUIM, ">$nome.$chr.$ini-$fim.PARTE_SNPS.M$map.SemCIGAR.tsv") or die "Não abriu o arquivo com parte dos SNPs! \n";

    foreach my $b_string(@haplo_ruim) {

	print OUTRUIM $b_string, "\n";

    } #foreach my $b_string(@haplo_ruim)

    close (OUTRUIM);
    #print "$cont_ruim\n\n"

} #while (my $pos = <POS>)

close (POS);
system ("mkdir $nome.PARTE_SNPS.M$map.SemCIGAR $nome.TODOS_SNPS.M$map.SemCIGAR $nome.CIGAR_RUIM.M$map.SemCIGAR");
system ("mv $nome.*.PARTE_SNPS.M$map.SemCIGAR.tsv $nome.PARTE_SNPS.M$map.SemCIGAR");
system ("mv $nome.*.TODOS_SNPS.M$map.SemCIGAR.tsv $nome.TODOS_SNPS.M$map.SemCIGAR");
system ("mv $nome.*.CIGAR_RUIM.M$map.SemCIGAR.tsv $nome.CIGAR_RUIM.M$map.SemCIGAR");
