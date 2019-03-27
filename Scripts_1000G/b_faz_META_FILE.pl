#!/usr/bin/perl -w

## Autora: Jaqueline Wang
## Mestranda do programa Interunidades em Bioinformática da USP

## Esse script serve para criar os arquivos de saída para cada microhaplotipo e cada populacao possivel

use strict;
use Getopt::Long;

my $help = 0;
my $num = 0;

GetOptions("help|h" => \$help,
	   "n=s" => \$num,
) or die "Erro ao pegar as opções! \n";

if ($help || !($num)) {die "\
$0: \
    Esse script recebe três arquivos, um é a lista de SNPs que compõem o microhaplótipo, o segundo é o vcf do 1000 Genomes e o terceiro é o número do micro-haplótipo.\
\
Parâmetros:\
     -h ou --help : Mostra as opções.\
     -n : numero do micro-haplótipo.\

\n";
}

my %haplos;
my $title = "Haplotipo\tTODOS";
my $total = "Frequencia";

open (POP, "/home/jaque/Desktop/Scripts_Martin/Arquivos/lista_pop.txt") or die "Não foi possível abrir o arquivo de lista de populacoes! \n";

open (OUT, ">M$num.meta_file.txt") or die "Não foi possível abrir o arquivo de saída! \n";

#Eliminamos a linha do TODOS!
<POP>;

#Abrimos o arquivo TODOS porque ele contém todos os haplótipos
open (TODOS, "Freq_haplotipo/M$num/M$num.TODOS.freq.txt") or die "Não foi possível abrir o arquivo TODOS! \n";


my $head = <TODOS>;
my $tot1 = 0;

while (my $line = <TODOS>) {	    
    chomp ($line);
    my ($hap, $freq) = split(/\t/, $line);
    
    $haplos{$hap} = $freq;
    $tot1 = $tot1 + $freq;
}

$total = $total."\t".$tot1;
	
close (TODOS);

while (my $line1 = <POP>) {
    chomp ($line1);
    my ($sigla, $pop) = split(/\t/, $line1);

    my %hash;
    $title = $title."\t".$pop;
 
    open (IN, "Freq_haplotipo/M$num/M$num.$pop.freq.txt") or die "Não foi possível abrir o arquivo $pop! \n";
    <IN>;
    my $tot2 = 0;
	
    while (my $line2 = <IN>) {
	chomp ($line2);
	my ($hap, $freq) = split(/\t/, $line2);
	
	$hash{$hap} = $freq;
	$tot2 = $tot2 + $freq;
    }
	
    close (IN);

    $total = $total."\t".$tot2;

    foreach my $key(keys(%haplos)) {
	
	if (exists($hash{$key})) {
	    $haplos{$key} = $haplos{$key}."\t".$hash{$key};

	}

	elsif (!exists($hash{$key})) {
	    $haplos{$key} = $haplos{$key}."\t0";

	}

    }

}

close (POP);
#close (OUT);

my $cont = 1;

print OUT "ID\t$title\n";
foreach my $key1(sort keys %haplos) {
    
    my $id;
    
    if ((length $cont) == 1) {
	$id = "M$num"."H0$cont";
    }

    else {
	$id = "M$num"."H$cont";
    }
    
    print OUT "$id\t$key1\t$haplos{$key1}\n";
    $cont += 1;
}
print OUT " \t$total\n";


