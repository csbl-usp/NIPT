#!/usr/bin/perl -w

#Autora: Jaqueline Wang
#Mestre do Programa Interunidades em Bioinformática - USP

#Script para extrair os haplótipos dos reads que tiveram um mapeamento bom e que a cobertura esteja contemplando todos os SNPs do micro-haplótipo.

#Devem ser passados 5 parâmetros:
#1 - Arquivo BAM
#2 - Qualidade do mapeamento
#3 - COM ou SEM CIGAR
#4 - Qualidade das BASES
#5 - Porcentagem de BASES cobertas

#Endereços importantes (MARTIN)
#Micro-haplótipos: /home/jaque/Desktop/Scripts_Martin/Arquivos/microhaplotipos.txt
#Arquivo BED: /home/jaque/Desktop/Scripts_Martin/Arquivos/SNPs.bed 
#Phred file: /home/jaque/Dektop/Scripts_Martin/Arquivos/phred.txt

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
Esse script recebe cinco entradas, o arquivo BAM, a qualidade de mapeamento, com ou sem CIGAR, qualidade das bases e porcentagem de bases cobertas. A saída são dois arquivos: \
MAIS_XX - contém os haplótipos com mais de XX% de bases cobertas. \
MENOS_XX - contém os haplótipos com mais de XX% de bases cobertas. \ 
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


my %qualidade;
open (QUAL, "/home/jaque/Desktop/Scripts_Martin/Arquivos/phred.txt") or die "Não foi possível abrir o arquivo com as qualidades das bases! \n";

while (my $line1 = <QUAL>) {
    chomp ($line1);
    my ($qual1, $qual2) = split (/\t/, $line1);

    if (!exists $qualidade{$qual1}) {
	$qualidade{$qual1} = $qual2;

    }
 
} #while (my $line1 = <QUAL>)

close (QUAL);


open (POS, "/home/jaque/Desktop/Scripts_Martin/Arquivos/microhaplotipos.txt") or die "Não foi possível obter os micro-haplótipos! \n";


while (my $pos = <POS>) {

    chomp ($pos);
    
    #Usamos esse match para armazenar o cromossomo e a posição inicial e final do micro-haplótipo
    
    $pos =~ m/(chr[0-9]+):(.*)-(.*)/;
    my $chr = $1;
    my $ini = $2;
    my $fim = $3;

    #Abrimos o arquivo com TODOS os SNPs cobertos
    open (INFILE, "$nome.TODOS_SNPS.M$map.$cigar/$nome.$chr.$ini-$fim.TODOS_SNPS.M$map.$cigar.tsv") or die "Não foi possível abrir o arquivo com os haplótipos! \n";
    
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
    open (BED, "/home/jaque/Desktop/Scripts_Martin/Arquivos/SNPs.bed") or die "Não abriu o BED! \n";
    
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
    open (INFILE1, "$nome.PARTE_SNPS.M$map.$cigar/$nome.$chr.$ini-$fim.PARTE_SNPS.M$map.$cigar.tsv") or die "Não foi possível abrir o arquivo com os haplótipos! \n";   

    
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

    #Arquivo com MAIS
    open (OUTMAIS, ">$nome.$chr.$ini-$fim.MAIS_$por.Q$score.M$map.$cigar.tsv") or die "Não foi possível abrir o arquivo de saída1! \n";

    #Arquivo com MENOS
    open (OUTMENOS, ">$nome.$chr.$ini-$fim.MENOS_$por.Q$score.M$map.$cigar.tsv") or die "Não foi possível abrir o arquivo de saída1! \n";    

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

system ("mkdir $nome.MAIS_$por.Q$score.M$map.$cigar $nome.MENOS_$por.Q$score.M$map.$cigar");
system ("mv $nome.*.MAIS_$por.Q$score.M$map.$cigar.tsv $nome.MAIS_$por.Q$score.M$map.$cigar");
system ("mv $nome.*.MENOS_$por.Q$score.M$map.$cigar.tsv $nome.MENOS_$por.Q$score.M$map.$cigar");
