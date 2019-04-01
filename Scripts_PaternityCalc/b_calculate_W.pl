#!/usr/bin/perl -w

#Author: Jaqueline Wang
#MsC in Bioinformatics Graduate Program - USP

#SCRIPT TO CALCULATE THE PROBABILITY OF PATERNITY USING THE FILES FROM ALLEGED FATHER, MOTHER AND PLASMA.

#The script receives 8 parameters
#1 - Mapping quality
#2 - YesCIGAR or NotCIGAR
#3 - Bases quality
#4 - Percentage of covered bases
#5- Coverage of alleged father and mother regions
#6 - Coverage of plasma regions
#7 - Population from 1000 Genomes
#8 - Trio number

#Inbalance parameters
#1 - Superior
#2 - Errors
#3 - Doubt between three

####################################################################################

use strict;
use Getopt::Long;

my $help = 0;
my $cobPL = 1000;
my $map = 20;
my $cigar = "NotCIGAR";
my $qual = 20;
my $por = 70;
my $cob = 20;
my $erros = 10;
my $superior = 80;
my $duvida = 35;
my $Genome = "ALL";
my $trio;

GetOptions("help|h" => \$help,
	   "m=s" => \$map,
	   "l=s" => \$cigar,
	   "q=s" => \$qual,
	   "p=s" => \$por,
	   "c=s" => \$cob,
	   "f=s" => \$cobPL,
	   "e=s" => \$erros,
	   "s=s" => \$superior,
	   "d=s" => \$duvida,
	   "g=s" => \$Genome,
	   "t=s" => \$trio
    ) or die "Failed to take the options! \n";

if ($help || !($trio)) {die "\
This script requires one input, the trio number. \
Other parameters have default values, but can be changed. \
The outputs is the probability of paternity and a report. \
\
Required parameters:\
	-t	Trio number \
\
Other parameters: \
	-h	Show the options \
	-m	Mapping quality of reads (default = 20) \
	-l	Uses or not the CIGAR info. If NOT (NotCIGAR), consider only (mis)matches. Other option is YesCIGAR (default = NotCIGAR) \
	-q	Bases quality (default = 20) \
	-p	Percentage of covered bases (default = 70) \
	-c	Coverage of alleged father and mother (default = 20) \
	-f	Coverage of plasma (default = 1000) \
	-g	Population from 1000 Genomes (default = ALL) \
\
Inbalance parameters: \
	-e	Limit for sequencing errors (default = 10) \
	-s	Limit to consider HOMOZYGOUS (default = 80) \
	-d	Limit to consider HETEROZYGOUS when there are more than 3 possibilities (default = 35) \
\n";
}

####################################################################################

my $nomeAF = "Trio$trio"."_SampleAF";
my $nomeM = "Trio$trio"."_SampleM";
my $nomeP = "Trio$trio"."_SampleP";

my $IPC = 1;
my $conta_m = 0;
my $tem_ff = 0;
my $tem_mut = 0;

####################################################################################

open (POS, "Files/Microhaplotypes.txt") or die "Failed to obtain the microhaplotypes! \n";

my $M = 1;

while (my $pos = <POS>) {
    chomp ($pos);
    print "$pos\n";

    #Usamos esse match para armazenar o cromossomo e a posição inicial e final do micro-haplótipo
    $pos =~ m/(chr[0-9]+):(.*)-(.*)/;
    my $chr = $1;
    my $ini = $2;
    my $fim = $3;
    
    ###################################################################################
    #Abrimos os arquivos de qualidade dos haplotipos e armazenamos as informacoes de cobertura e porcentagem de cada haplotipo encontrado
    open (INFILE, "BAM/Trio$trio/$nomeP/$nomeP.haplotipos.MORE_$por.Q$qual.M$map.$cigar.tsv") or die "Failed to open HAPLOTYPES file! \n";
    
    my %cobertura;
    my %porcentagem;
    my $n_reads = 0;
    
    while (my $line1 = <INFILE>) {
	chomp ($line1);
	next unless $line1 =~ m/^pos/;
	next unless $line1 !~ m/Discard/;
	my @dados = split(/\t/, $line1);

	if ($pos eq $dados[0]) {
	    $porcentagem{$dados[1]} = $dados[3];
	    $cobertura{$dados[1]} = $dados[2];
	    $n_reads = $n_reads + $cobertura{$dados[1]};
	}
    }
    close (INFILE);

    ####################################################################################
    #Buscamos os haplotipos do pai no arquivo do PAI
    open (PAI, "Genotypes_$nomeAF.M$map.$cigar.Q$qual.P$por.C$cob.E$erros.S$superior.D$duvida.txt") or die "Failed to open the AF GENOTYPES file!\n";
    $/ = "\n\n";
    my %haplo_pai;
    
    while (my $line2 = <PAI>) {
	chomp ($line2);

	if ($line2 =~ m/$pos/g) {
	    my @array = split(/\n/, $line2);
	    
	    foreach my $key(@array) {
		if ($key =~ m/^[ATCG]+/g) {
		    my ($haplo, $porce, $cober) = split(/\t/, $key);
		    $haplo_pai{$haplo} = $porce."\t".$cober;
		}
	    }
	}
    }
    close (PAI);

    ####################################################################################
    #Buscamos os haplotipos da mae no arquivo da MAE
    open (MAE, "Genotypes_$nomeM.M$map.$cigar.Q$qual.P$por.C$cob.E$erros.S$superior.D$duvida.txt") or die "Failed to open the M GENOTYPES file! \n";
    my %haplo_mae;
    
    while (my $line3 = <MAE>) {
	chomp ($line3);
	
	if ($line3 =~ m/$pos/g) {
	    my @array = split(/\n/, $line3);
	    
	    foreach my $key(@array) {

		if ($key =~ m/^[ATCG]+/g) {
		    my ($haplo, $porce, $cober) = split(/\t/, $key);
		    $haplo_mae{$haplo} = $porce."\t".$cober;
		}
	    }
	}
    }
    close (MAE);

    $/ = "\n";
    
    #Buscamos os haplotipos do feto que podem ter sido herdados do pai
    my %herda_pai;
    my %haplo_erro;
    my %mae_plasma;
    my %problemas;
    
    foreach my $key1(keys(%porcentagem)) {

	if (exists($haplo_mae{$key1})) {
	    if (($porcentagem{$key1} > 12) && ($n_reads >= $cob)) {
		my $num = sprintf "%.2f", $porcentagem{$key1};
		$mae_plasma{$key1} = "$num%\t$cobertura{$key1}X";
	    }
	}

	elsif(!exists($haplo_mae{$key1})) {
	    if (($porcentagem{$key1} >= 1) && ($porcentagem{$key1} <= 12) && ($n_reads >= $cob)) {
		
		if (exists($haplo_pai{$key1})) {
		    my $num = sprintf "%.2f", $porcentagem{$key1};
		    $herda_pai{$key1} = "$num%\t$cobertura{$key1}X";
		}

		elsif (!exists($haplo_pai{$key1})) {
		    my $num = sprintf "%.2f", $porcentagem{$key1};
		    $haplo_erro{$key1} = "$num%\t$cobertura{$key1}X";
		}
	    }

	    if (($porcentagem{$key1} > 12) && ($n_reads >= $cob)) {
		my $num = sprintf "%.2f", $porcentagem{$key1};
		$problemas{$key1} = "$num%\t$cobertura{$key1}X";

	    }
	}
    }
	
    #Verificamos se existem haplotipos comuns entre o suposto pai e a mae
    if (((scalar keys %haplo_pai) > 0) && ((scalar keys %haplo_mae) > 0) && ($n_reads>= $cob) && ((scalar keys %haplo_mae) == (scalar keys %mae_plasma)) && ((scalar keys %problemas) == 0)) {
	
	print "HAPLO ALLEGED FATHER!\n";
	foreach my $keyp(keys(%haplo_pai)) {
	    print "$keyp\t$haplo_pai{$keyp}\n";
	}

	print "\n";
	print "HAPLO MOTHER!\n";

	foreach my $keym(keys(%haplo_mae)) {
	    print "$keym\t$haplo_mae{$keym}\n";
	}

	print "HAPLO MOTHER-PLASMA\n";

	foreach my $keymp(keys(%mae_plasma)) {
	    print "$keymp\t$mae_plasma{$keymp}\n";
	}

	print "\n";
	print "INHERIT ALLEGED FATHER! \n";

	foreach my $keyf(keys(%herda_pai)) {
	    print "$keyf\t$herda_pai{$keyf}\n";
	}

	print "\n";
	print "SEQ ERROR! \n";
	my $dist_erro = 0;
	
	foreach my $keye(keys(%haplo_erro)) {
	    my @e_base = $keye =~ m/[ATCG]/g;
	    my $dist = scalar(@e_base);

	    foreach my $key1(keys(%herda_pai)) {
		my $cont1 = 0;
		my $cont2 = 0;
		my @p_base = $key1 =~ m/[ATCG]/g;
				
		while ($cont1 < scalar(@p_base)) {
		    if ($p_base[$cont1] ne $e_base[$cont1]) {
			$cont1 += 1;
			$cont2 += 1;
		    }
		    
		    else {
			$cont1 += 1;
		    }
		}
		
		if ($dist > $cont2) {
		    $dist = $cont2;
		}
	    }
	    
	    foreach my $key2(keys(%haplo_mae)) {
		my $cont1 = 0;
		my $cont2 = 0;
		my @m_base = $key2 =~ m/[ATCG]/g;
			
		while ($cont1 < scalar(@m_base)) {
		    if ($m_base[$cont1] ne $e_base[$cont1]) {
			$cont1 += 1;
			$cont2 += 1;
		    }
		    
		    else {
			$cont1 += 1;
		    }
		}
		
		if ($dist > $cont2) {
		    $dist = $cont2;
		}
	    }
    
	    print "$keye\t$haplo_erro{$keye}\t$dist\n";
	    
	    if ($dist_erro < $dist) {
		$dist_erro = $dist;
	    }
	}
	print "\n";

	if ((scalar keys %herda_pai) > 0) {
	    $tem_ff += 1;
	}

	my $ev1 = 0;
	foreach my $key1(keys(%haplo_pai)) {
	    if (exists($haplo_mae{$key1})) {
		$ev1 += 1;
	    }
	}
	
	my $ev = 0;  
	
	if ((scalar keys %herda_pai) == 0) {
	    if ($ev1 >= 1) {
		$ev += 0.5;
	    }

	    else {
		$tem_mut += 1;
	    }
	}
	
	elsif ((scalar keys %herda_pai) == 1) {
	    $ev += 1;
	}

	elsif ((scalar keys %herda_pai) >= 2) {
	    $ev += 1;
	}

	print "RESUME\n";
	print "EV = $ev\n";
	$conta_m += 1;
	
	##Calcular o IP!!
	my $n_h_mae = scalar keys %haplo_mae;
	my $n_h_pai = scalar keys %haplo_pai;
	my $n_h_feto = scalar keys %herda_pai;
	my $IP = 0;
	my $meta_file;
	
	if ($M < 10) {
	    $meta_file = "M0".$M.".meta_file.txt";
	}

	elsif ($M >= 10) {
	    $meta_file = "M".$M.".meta_file.txt";
	}

	open (META, "Haplotypes/$meta_file") or die "Failed to open meta file! \n";
	my $meta_head = <META>;
	chomp ($meta_head);
	my @dado_head = split(/\t/, $meta_head);
	my $n_col = 0;
	my $cont_g = 0;
	
	foreach my $key_m(@dado_head) {
	    if ($key_m eq $Genome) {
		$n_col = $cont_g;
	    }
	    $cont_g += 1;
	}
	
	my %meta_haplo;
	my $total = 0;
	
	while (my $line4 = <META>) {
	    if ($line4 =~ m/M[0-9]{2}H[0-9]{2}/g) {
		chomp ($line4);
		my @dado_meta = split(/\t/, $line4);
		$meta_haplo{$dado_meta[1]} = $dado_meta[$n_col];
	    }

	    elsif ($line4 =~ m/\tFrequency\t/g) {
		chomp($line4);
		my @dado_meta = split(/\t/, $line4);
		$total += $dado_meta[$n_col];
	    }
	}
	close (META);

	if ($ev == 0) {
	    print "Case MUTATION\n";
	    $IP += 10**(-8);    
	}

	elsif ($ev > 0) {
	    my $count1 = 0;
	    my @ab;
	    my $haplo_b;
  
	    foreach my $key1(keys(%haplo_mae)) {
		if (exists($haplo_pai{$key1})) {
		    $count1 += 1;
		}

		elsif(!exists($haplo_pai{$key1})) {
		    $haplo_b = $key1;
		}
		push @ab, $key1;
	    } #foreach my $key1(keys(%haplo_mae)

	    if ($n_h_feto == 0) {
		my $a = 0;
		my $b = 0;
		    
		if (($n_h_pai == 2) && ($n_h_mae == 2)) {
		    if ($count1 == 1) {

			if ($ab[0] ne $haplo_b) {
			    $a += $meta_haplo{$ab[0]}/$total;
			    $b += $meta_haplo{$ab[1]}/$total;
			}

			elsif ($ab[1] ne $haplo_b) {
			    $a = $meta_haplo{$ab[1]}/$total;
			    $b = $meta_haplo{$ab[0]}/$total;		
			}

			if ($a == 0) {
			    $a += 1/($total+1);
			}

			if ($b == 0) {
			    $b += 1/($total+1);
			}
			
			#print "A = $a\n";
			#print "B = $b\n";
			
			print "Case M=AB, C=AB, AF=AC ou M=AB, C=A, AF=AC\n";

			my $IP1 = (1/(2*$a));
			my $IP2 = (1/(2*($a + $b)));

			if ($IP1 <= $IP2) {
			    $IP += $IP1;

			}

			elsif ($IP1 > $IP2) {
			    $IP += $IP2;

			}
		    }

		    elsif ($count1 == 2) {

			$a += $meta_haplo{$ab[0]}/$total;
			$b += $meta_haplo{$ab[1]}/$total;

			if ($a == 0) {
			    $a += 1/($total+1);
			}

			if ($b == 0) {
			    $b += 1/($total+1);
			}
			
			#print "A = $a\n";
			#print "B = $b\n";
			
			print "Case M=AB, C=AB, AF=AB ou M=AB, C=A, AF=AB\n";

			if ($a < $b) {
			    $IP += (1/(2*$b));
			}

			elsif ($a == $b) {
			    $IP += (1/($a + $b));
			}

			elsif ($a > $b) {
			    $IP += (1/(2*$a));
			}
		    }
		}

		elsif (($n_h_pai == 1) && ($n_h_mae == 2)) {
   
		    if ($ab[0] ne $haplo_b) {
			$a += $meta_haplo{$ab[0]}/$total;
			$b += $meta_haplo{$ab[1]}/$total;
		    }

		    elsif ($ab[1] ne $haplo_b) {
			$a = $meta_haplo{$ab[1]}/$total;
			$b = $meta_haplo{$ab[0]}/$total;		
		    }

		    if ($a == 0) {
			$a += 1/($total+1);
		    }

		    if ($b == 0) {
			$b += 1/($total+1);
		    }
		    
		    #print "A = $a\n";
		    #print "B = $b\n";
		    
		    if ($count1 == 1) {
			print "Case M=AB, C=AB, AF=A ou M=AB, C=A, AF=A\n";

			my $IP1 = (1/($a));
			my $IP2 = (1/($a + $b));

			if ($IP1 <= $IP2) {
			    $IP += $IP1;
			}

			elsif ($IP1 > $IP2) {
			    $IP += $IP2;
			}
		    }
		}

		elsif (($n_h_pai == 2) && ($n_h_mae == 1)) {

		    $a += $meta_haplo{$ab[0]}/$total;

		    if ($a == 0) {
			$a += 1/($total+1);
		    }
		    
		    #print "A = $a\n";
		    
		    if ($count1 == 1) {
			print "Case M=A, C=A, AF=AB\n";
			$IP += (1/(2*$a));
		    }		    
		}

		elsif (($n_h_pai == 1) && ($n_h_mae == 1)) {

		    $a += $meta_haplo{$ab[0]}/$total;

		    if ($a == 0) {
			$a += 1/($total+1);
		    }
		    
		    #print "A = $a\n";
		    
		    if ($count1 == 1) {
			print "Case M=A, C=A, AF=A\n";
			$IP += (1/$a);
		    }
		}
		
	    } #if ($n_h_feto == 0)

	    elsif ($n_h_feto == 1) {

		my $a = 0;
		
		foreach my $key2(keys(%herda_pai)) {
		    $a += $meta_haplo{$key2}/$total;
		}

		if ($a == 0) {
		    $a += 1/($total+1);
		}
		
		#print "A = $a\n";
		if (($n_h_pai == 2) && ($n_h_mae == 2)) {
		    if ($count1 == 0) {
			print "Case M=BD, C=AB, SP=AC\n";
			$IP += (1/(2*$a));
		    }

		    elsif ($count1 == 1) {
			print "Case M=BC, C=AB, AF=AC ou M=BC, C=AB, AF=AB\n";
			$IP += (1/(2*$a));
		    }
		}

		elsif (($n_h_pai == 1) && ($n_h_mae == 2)) {
		    if ($count1 == 0) {
			print "Case M=BC, C=AB, AF=A\n";
			$IP += (1/$a);
		    }
		}

		elsif (($n_h_pai == 2) && ($n_h_mae == 1)) {
		    if ($count1 == 0) {
			print "Case M=B, C=AB, AF=AC\n";
			$IP += (1/(2*$a));
		    }

		    elsif ($count1 == 1) {
			print "Case M=B, C=AB, AF=AB\n";
			$IP += (1/(2*$a));
		    }
		}

		elsif (($n_h_pai == 1) && ($n_h_mae == 1)) {
		    if ($count1 == 0) {
			print "Case M=B, C=AB, AF=A\n";
			$IP += (1/$a);
		    }
		}
	    } #elsif ($n_h_feto == 1)

	    elsif ($n_h_feto == 2) {

		my $a1 = 0;
		my $a = 0;
		    
		foreach my $key2(keys(%herda_pai)) {
		    $a1 += $meta_haplo{$key2}/$total;
		
		    if ($a1 == 0) {
			$a1 += 1/($total+1);
		    }

		    if ($a1 > $a) {
			$a = $a1;
		    }
		}
		
		#print "A = $a\n";
		if (($n_h_pai == 2) && ($n_h_mae == 2)) {
		    if ($count1 == 0) {
			print "Case M=BD, C=AB, AF=AC\n";
			$IP += (1/(2*$a));
		    }

		    elsif ($count1 == 1) {
			print "Case M=BC, C=AB, AF=AC ou M=BC, C=AB, AF=AB\n";
			$IP += (1/(2*$a));
		    }
		}

		elsif (($n_h_pai == 1) && ($n_h_mae == 2)) {
		    if ($count1 == 0) {
			print "Case M=BC, C=AB, AF=A\n";
			$IP += (1/$a);
		    }
		}

		elsif (($n_h_pai == 2) && ($n_h_mae == 1)) {
		    if ($count1 == 0) {
			print "Case M=B, C=AB, AF=AC\n";
			$IP += (1/(2*$a));
		    }

		    elsif ($count1 == 1) {
			print "Case M=B, C=AB, AF=AB\n";
			$IP += (1/(2*$a));
		    }
		}

		elsif (($n_h_pai == 1) && ($n_h_mae == 1)) {
		    if ($count1 == 0) {
			print "Case M=B, C=AB, AF=A\n";
			$IP += (1/$a);
		    }
		}
	    } #elsif ($n_h_feto == 2)
	}
	
	$IPC = $IPC*$IP;
	print "IP = $IP\n";
    }
    print "\n";
    $M += 1;
} #while (my $pos = <POS>)

close (POS);

my $W = $IPC/($IPC + 1);
print "REPORT\n";
print "Good microhaplotypes = $conta_m\n";
print "Microhaplotypes with FF = $tem_ff\n";
print "Microhaplotypes with mutation = $tem_mut\n";
print "IPC = $IPC\n";
print "W = $W\n";
print "\n";

print "QUALITY INFORMATION\n";
print "Mapping quality of reads = Q$map\n";
print "CIGAR string = $cigar\n";
print "Bases quality = Q$qual\n";
print "Percentage of covered bases = $por%\n";
print "Coverage of ALLEGED FATHER and MOTHER = $cob\n";
print "Coverage of PLASMA = $cobPL\n";
print "Population from 1000 Genomes = $Genome\n";
print "Limit for sequencing errors = $erros%\n";
print "Limit to consider HOMOZYGOUS = $superior%\n";
print "Limit to consider HETEROZYGOUS when ther are more than 3 possibilities = $duvida%\n";
print "\n";
