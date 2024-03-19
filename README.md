# MBMDRc_R
 Implementation in R of MBMDRc


This library provides a clean and easy implementation of MBMDRc ( https://pubmed.ncbi.nlm.nih.gov/33602124/ ) in the case of continuous phenotype. 
There are two different scripts for MBMDR: 1d (only main effect) and 2d (only interaction). 

After installation, it is possible to run the code with the appropriate parameters: 

* p_threshold, the p-value threshold: all the SNP (or SNP-pairs if 2d) under the p-value are considered 
* is_0_consider: {0;1}, if 0, the cells of the HLO matrix marked with 0 are also considered in the calculation; otherwise, only H and L. 
* data_path: a string for the data path
* result_path: a string for where you want to save the results
* model: the name and location (in data path if not specified) of the MBMDR model
* output: the name and location (in data path if not specified) of the hits found my MBMDR
* data: the name and location of the input MBMDR data

Example with MockData: 

``` ./MBMDRc_1d_script.R 0.05 0 ../Data/1d/ ../../Result/1d/ mbmdr_model_1d.txt mbmdr_output_1d.txt MockData_1d.txt ```

or, for 2d 

``` ./MBMDRc_2d_script.R 0.05 0  ../Data/2d/ ../../Result/2d/ MBMDRinput_models.txt mbmdr_output.txt MockData.txt ```

Alternatively, it is possible to open the .R file in Rstudio and manually set the parameters to their value.

The required inputs have the following shape: 

* mbmdr_model_1d, (MBMDRinput_models if 2d) with the subject, average continuous trait and HL0 matrix for each hit
* mbmdr_output_1d, (mbmdr_output if 2d), with the top hits (either SNp or SNP pairs), togheter with the chi square statistic and the p-value
* MockData: data simulating a tabular matrix with individuals on the rows and SNPs on the column, with the phenotype as the first column.

For the wrapper (Wrapped_MBMDRc.R), the parameters are the same as before for 1d and 2d with only the addition of the parameter d_level, as last parameter, that can take the value "1" (MBMDRc-1d) or "2" (MBMDRc-2d). A typical call to the wrapper look a bit something like this: 

``` Wrapper_MBMDRc("0.05", "0", "../Data/2d/", "../../Result/2d/", "MBMDRinput_models.txt", "mbmdr_output.txt", "MockData.txt", "2") ```

