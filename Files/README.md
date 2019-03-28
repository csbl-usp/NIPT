# Files

The files located in this folder are used to support the analysis performed by the scripts.

**1. Microhaplotypes.txt** - This file is used to keep the microhaplotype genomic position, you have to put the chromosome number, start position and end position in the following format:
  ```
  chr4:7447228-7447353
  ```
  
**2. SNPs.bed** - This file is used to keep the positions of analyzed SNPs, and also their ID (rs number). It has to be in 7 columns bed format.
  ```
  chr4	7447228	7447228	rs11721645	chr4	7447228	7447228
  ```
  
**3. integrated_call_samples_v3.20130502.ALL.panel** - This file was obtained from 1000 Genomes, and it keeps the information of each sample of the database. It is used to obtain the population, the super population and the gender of each sample.

**4. phred.txt** - This file is used to obtain the character that represents each quality value in phred scale.

**5. pop_list.txt** - This file is used to obtain the label of each super-population and population.
