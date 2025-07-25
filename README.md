# WHO_NI_Tutorial

This tutorial provides guidance on making non-inferiority assessments from experimental hut trial (EHT) data. Code is provided in both R and STATA. You should begin by downloading the files, unzipping the folder, and moving the folder to your desired location on your computer. To use the R code, you should first download [R](https://cran.r-project.org/) and [Rstudio](https://www.rstudio.com). We recommend opening the R Project file 'WHO_NI_tutorial.RProj' first: this means that your R session in RStudio should begin with the working directory set to the folder containing the relevant files. 

This work supports the WHO Report [Technical consultation to assess comparative efficacy of vector control products: meeting report, 5 and 9 June 2023](https://www.who.int/publications/i/item/9789240078659).

The R script 'power_calculator_user_script.R' contains a simulation-based method for estimating study power. The script contains a lengthy explanation of how a study should be characterised, **which should be read carefully**. The methods here have been adapted from [earlier work](https://github.com/JDChallenger/EHT_Visualise). This work built on a study by [Johnston et al.](https://doi.org/10.1111/2041-210X.12306). The main difference between this work and the earlier work is that the simulated datasets are being analysed with generalised linear regression models that contain only fixed effects (no random effects). This is consistent with the tutorial hosted on this repository. Additionally, the methodology used here allows for a non-inferiority margin that varies with the performance of the active comparator product.

Before power can be calculated, the design of the EHT must be specified. As outlined [here](https://doi.org/10.1016/j.crpvbd.2023.100115), the number of mosquitoes entering the hut and the level of variability present in the assay will both influence power (as will the trial duration).  The script 'power_calculator_user_script.R' provides details on how the scenario in question should be specified. The functions that this script uses are contained in 'power_calculator_functions_FE.R'. 

## Update (February 2025): Community Studies
This repository has been updated, to also include a tutorial on assessing field-aged nets, collected during community studies. All the materials pertaining to community studies (R code and the PDF version of the tutorial) can be found in the folder 'Community_Studies'. A hyperlink to the relevant WHO documentation will appear here when it is available

## Update (May 2025): Regeneration Studies for ITNs
This repository has been updated to include guidance on estimating the regeneration time for ITNs, after washing. The relevant R code can be found in the folder 'Regeneration_Time'. Here, we include code for analysing data collected from a bioassay (measuring mosquito mortality) and from a chemical analysis. The Excel sheets in this folder contain some simulated datasets, to allow users to practise carrying out the procedure, before applying it to their own dataset(s). Please note: these simulated datasets were updated in July 2025, to provide greater realism.

The relevant WHO documentation for these studies can be found [here](https://extranet.who.int/prequal/node/29912).
