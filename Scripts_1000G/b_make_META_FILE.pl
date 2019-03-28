#!/usr/bin/perl -w

#Author: Jaqueline Wang
#MsC in Bioinformatics Graduate Program - USP

#SCRIPT TO CREATE THE META FILE WITH ALL THE POPULATIONS

use strict;
use Getopt::Long;

my $help = 0;
my $num = 0;

GetOptions("help|h" => \$help,
	   "n=s" => \$num,
) or die "Failed to take the options! \n";

if ($help || !($num)) {die "\
This script receives one parameter. \
\
Parameter: \
     -h	Show the options \
     -n	Microhaplotype number \
\n";
}

my %haplos;
my $title = "Haplotype\tALL";
my $total = "Frequency";

open (POP, "Files/pop_list.txt") or die "Failed to open POP_LIST file! \n";

open (OUT, ">M$num.meta_file.txt") or die "Failed to open META FILE! \n";

#Eliminamos a linha do TODOS!
<POP>;

#Abrimos o arquivo TODOS porque ele contém todos os haplótipos
open (TODOS, "1000G_test/Freq_haplo/M$num/M$num.ALL.freq.txt") or die "Failed o open ALL file! \n";


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
 
    open (IN, "1000G_test/Freq_haplo/M$num/M$num.$pop.freq.txt") or die "Failed to open $pop file! \n";
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

close (OUT);
