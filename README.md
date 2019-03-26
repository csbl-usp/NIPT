# NIPT
SCRIPTS PARA CALCULAR PROBABILIDADE DE PATERNIDADE

1 - Para fazer o cálculo da probabilidade de paternidade:

a - Entrar na pasta onde estão os arquivos BAM do Suposto Pai, Mãe e Plasma.
b - Chamar o script Paternity_Calc.pl, linha de comando - /home/martin/Desktop/Scripts_Martin/Paternity_Calc.pl -X BAM_suposto_pai -Y BAM_mae -Z BAM_plasma
c - Caso queira, pode mudar os parâmetros.

2 - Para incluir um novo micro-haplótipo nas análises:

a - Entrar na pasta /home/martin/Desktop/Scripts_Martin/teste_1000G
b - Colocar o vcf do novo micro-haplótipo na pasta VCF. (Ex de nome do vcf: 4.7447228-7447353.ALL.chr4.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf)
c - Colocar a lista de SNPs do novo micro-haplótipo na pasta MICRO. (Ex de nome da lista de SNPs: lista4.7447228-7447353)
d - Chamar o script 1000G.pl, linha de comando - /home/martin/Desktop/Scripts_Martin/1000G.pl -v VCF/nome_vcf -m MICRO/nome_lista -n NUM_MICRO
e - A saída será um arquivo chamado MXX_meta_file.txt (onde XX = NUM_MICRO).
f - Colocar o arquivo MXX_meta_file.txt na pasta /home/martin/Desktop/Scripts_Martin/Haplotipos
g - Abrir o arquivo /home/martin/Desktop/Scripts_Martin/Arquivos/microhaplotipos.txt, e adicionar uma linha com o chromossomo, posição inicial e posição final do micro-haplótipo. Ex: chr4:7447228-7447353
h - Adicionar no arquivo /home/martin/Desktop/Scripts_Martin/SNPs.bed os SNPs que fazem parte do novo micro-haplótipo no formato BED de 7 colunas.
