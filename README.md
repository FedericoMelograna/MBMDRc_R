# MBMDRc_R
 Implementation in R of MBMDRc


This library provides a clean and easy implementation of MBMDRc ( https://pubmed.ncbi.nlm.nih.gov/33602124/ ) in the case of continuous phenotype. 
There are two different scripts for MBMDR: 1d (only main effect) and 2d (only interaction). 

After installation, it is possible to run the code with the appropriate parameters: 

* p_threshold, the p-value threshold: all the SNP (or SNP-pairs if 2d) under the p-value are considered 
* is_0_consider: {0;1}, if 0, the cells of the HLO matrix marked with 0 are also considered in the calculation; otherwise, only H and L. 

Example: 

``` ./MBMDRc_1d_script.R 0.05 0 ```

or, for 2d 

``` ./MBMDRc_2d_script.R 0.05 0  ```

Alternatively, it is possible to open the .R file in Rstudio and manually set the parameters to their value.

The required inputs are: 

* mbmdr_model_1d, (MBMDRinput_models if 2d) with the subject, average continuous trait and HL0 matrix for each hit
* mbmdr_output_1d, (mbmdr_output if 2d), with the top hits (either SNp or SNP pairs), togheter with the chi square statistic and the p-value
* MockData: data simulating a tabular matrix with individuals on the rows and SNPs on the column, with the phenotype as the first column.

 
