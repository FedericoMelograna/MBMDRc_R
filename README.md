# MBMDRc_R
 Implementation in R of MBMDRc


This library provides a clean and easy implementation of MBMDRc ( https://pubmed.ncbi.nlm.nih.gov/33602124/ ) in the case of continuous phenotype. 
There are two different scripts for MBMDR: 1d (only main effect) and 2d (only interaction). 

After installation, it is possible to run the code with the appropriate parameters: 

* p_threshold, the p-value threshold: all the SNP (or SNP-pairs if 2d) under the p-value are considered 
* is_0_consider: {0;1}, if 0, the cells of the HLO matrix marked with 0 are also considered in the calculation; otherwise, only H and L. 
* model: the name and location of the MBMDR model
* data: the name and location of the input MBMDR data
* chisq: the name and location of the hits found my MBMDR
* d_level: {1,2}, with 1 if we want to run MBMDR 1d and 2 if we want to run MBMDR 2d

Example with MockData: 

``` ./mbmdrc_wrapper.R 0.05 0 ../Data/1d/mbmdr_model_1d.txt ../Data/1d/MockData_1d.txt ../Data/1d/mbmdr_output_1d.txt 1 ```

or, for 2d 

``` ./mbmdrc_wrapper.R 0.05 0 ../Data/2d/MBMDRinput_models.txt ../Data/2d/MockData_1d.txt ../Data/2d/mbmdr_output.txt 2 ```


The input of the simulated example are as follows: 

* mbmdr_model_1d, (MBMDRinput_models if 2d) with the subject, average continuous trait and HL0 matrix for each hit
* mbmdr_output_1d, (mbmdr_output if 2d), with the top hits (either SNp or SNP pairs), togheter with the chi square statistic and the p-value
* MockData: data simulating a tabular matrix with individuals on the rows and SNPs on the column, with the phenotype as the first column.

Alternatively, it is possible to open the .R file in Rstudio and manually set the parameters to their value. A typical call to the wrapper looks like this: 

```  
Wrapper_MBMDRc(p_threshold = "0.05", is_0_considered = "0",  model = "../Data/1d/mbmdr_model_1d.txt",data = "../Data/1d/MockData_1d.txt",  chisq = "../Data/1d/mbmdr_output_1d.txt", d_level = "1")
Wrapper_MBMDRc(p_threshold = "0.05", is_0_considered = "0",  model = "../Data/2d/MBMDRinput_models.txt",data = "../Data/2d/MockData.txt",  chisq = "../Data/2d/mbmdr_output.txt", d_level = "2")
```

