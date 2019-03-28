# Non Invasive Prenatal Testing (NIPT)

### Scripts to calculate the probability of paternity using NGS data of mother, plasma and alleged father

**1. Download NIPT package**

  - Download the NIPT package.
  
  - Using the terminal, change the files permissions with the following command.
    ```
    chmod 755 -R NIPT_master
    ```
  
**2. Download the Samtools package**

  - Download the Samtools package, versions 1.3.1.
    ```
    https://sourceforge.net/projects/samtools/files/samtools/1.3.1/
    ```
  
  - Unzip the Samtools package and put inside the **Samtools** folder.
  
  - Inside the folder **Samtools/samtools-1.3.1**, run the following commands:
    ```
    ./configure
    make
    ```    

**2. To include a new microhaplotype for analisys**
  
  - Inside the folder **1000G_test**, there are two folders, onde named **VCF** and other named **MICRO**.
  
  - Put the vcf (obtained from 1000 Genomes) in the folder **VCF**. Example of vcf name: 
      ```
      4.7447228-7447353.ALL.chr4.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf
      ```
  
  - Put the microhaplotype's list of SNPs in the folder **MICRO**. Example of list name:
      ```
      list4.7447228-7447353
      ```
  
  - Call the script 1000G.pl (inside the NIPT folder).
    ```
    ./1000G.pl -v 1000G_test/VCF/vcf_name -m 1000G_test/MICRO/list_name -n microhaplotype_number
    ```
    Obs: The microhaplotype number has to be sequencial and composed of two numbers.
    
  - The output is a file named MXX_meta_file.txt (where XX = microhaplotype_number).
  
  - Put the file MXX_meta_file.txt in the folder **Haplotypes**.
  
  - There is a file named **Microhaplotypes.txt** inside the folder **Files**. Add a line in the file with chromosome, start position and end position (chr4:7447228-7447353).
  
  - There is a file named **SNPs.bed** inside the folder **Files**. Put the SNPs that compose the microhaplotype in bed format (7 columns).


**3. To calculate the probability of paternity**

  - Go to the directory where the mother, the plasma and the alleged father's BAM files are located.
  
  - Call the script Paternity_Calc.pl 
    ```
    Paternity_Calc.pl -X BAM_alleged_father -Y BAM_mother -Z BAM_plasma
    ```
  
  - We settled the default parameters, but it can be changed.


