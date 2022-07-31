
# Copula based Cox proportional hazards models for dependent censoring
<!-- output:
  pdf_document: default
  html_document: default
-->
<!--HOW TO COMPLETE THIS FORM:-->

<!--
1. Checkboxes in this document appear as follows: 

- [ ] This is a checkbox 

To check a checkbox, replace [ ] by [x], as follows: 

- [x] This is a checked checkbox 

Note that older versions of RStudio (versions lower than 1.3) may not create a formatted checkbox but will leave the original characters, i.e., literally "[ ]" or "[x]". It's fine to submit a PDF in this form.
 
2. For text answers, simply type the relevant text in the areas indicated. A blank line starts a new paragraph. 
 
3. Comments (like these instructions) provide additional instructions throughout the form. There is no need to remove them; they will not appear in the compiled document. 

4. If you are comfortable with Markdown syntax, you may choose to include any Markdown-compliant formatting in the form. For example, you may wish to include R code chunks and compile this document in R Markdown.
-->

This form documents the artifacts associated with the article (i.e., the data and code supporting the computational findings) and describes how to reproduce the findings.


# 1. Data

<!-- 
[ ] This paper does not involve analysis of external data (i.e., no data are used or the only data are generated by the authors via simulation in their code).
-->

<!--
If box above is checked and if no simulated/synthetic data files are provided by the authors, please skip directly to the Code section. Otherwise, continue.
-->

- [x] I certify that the authors of the manuscript have legitimate access to and permission to use the data used in this manuscript.

<!-- If data are simulated using random number generation, please be sure to set the random number seed in the code you provide -->

## Abstract

A follicular cell lymphoma data set given in the book of Pintilie (2006) is used to illustrate the performance of the proposed semiparametric copula model. The data given in the book consist of 541 patients with early disease stage (I or II) and treated with radiation alone (RT) or with radiation and chemotherapy (CMT). The endpoints of interest are what comes first: relapse of the disease or death in remission. The data are freely available in the $R$ package **randomForestSRC**. We also studied the performance of the proposed method using extensive simulations.


<!--

Provide a short (< 100 words), high-level description of the data
-->

## Availability


- [x] Data **are** publicly available.

<!-- [ ] Data **cannot be made** publicly available.

If the data are publicly available, see the *Publicly available data* section. Otherwise, see the *Non-publicly available data* section, below.-->

### Publicly available data

- [x] Data are available online at: https://cran.r-project.org/web/packages/randomForestSRC/




<!-- If data are available by request to the authors or some other data owner, please make sure to explain the process of requesting access to the data. 

### Non-publicly available data -->

<!--
The Journal of the American Statistical Association requires authors to make data accompanying their papers available to the scientific community except in cases where: 1) public sharing of data would be impossible, 2) suitable synthetic data are provided which allow the main analyses to be replicated (recognizing that results may differ from the "real" data analyses), and 3) the scientific value of the results and methods outweigh the lack of reproducibility.

Please discuss the lack of publicly available data. For example:
-	why data sharing is not possible,
-	what synthetic data are provided, and 
-	why the value of the paper's scientific contribution outweighs the lack of reproducibility.
-->

## Description

The data are available in the $R$ package **randomForestSRC** by data name **follic**. Hence, all the descriptions about the data can be obtained by typing **?follic** in $R$. 

<!-- 
OPTIONAL: Provide any additional details that would be helpful in understanding the data. If relevant, please provide unique identifier/DOI/version information and/or license/terms of use.
-->

# 2. Code

## Abstract
All of the data simulations, processing and analysis for this paper were done in $R$.  The code provided includes all the scripts necessary to replicate all the tables and figures in the paper, for both real and simulated data. 


<!--
Provide a short (< 100 words), high-level description of the code. If necessary, more details can be provided in files that accompany the code. If no code is provided, please state this and say why (e.g., if the paper contains no computational work).
-->

## Description

Simulations in the manuscript were conducted using R version 3.6.0. All of the R scripts used in the paper are available in a public repository on GitHub [url blinded].  The MIT license applies to all code, and no permissions are required to access the code. 

To run the code provided, the following $R$ libraries should be installed. The numbers in the parentheses indicate the version of the packages used when performing the analysis presented in the manuscript. **pbivnorm (0.6.0),  survival (3.2-13), copula (1.0-1)**.


<!--
### Code format(s)

-->
<!--
Check all that apply
- [ ] Script files
    - [ ] R
    - [ ] Python
    - [ ] Matlab
    - [ ] Other: 
- [ ] Package
    - [ ] R
    - [ ] Python
    - [ ] MATLAB toolbox
    - [ ] Other: 
- [ ] Reproducible report 
    - [ ] R Markdown
    - [ ] Jupyter notebook
    - [ ] Other:
- [ ] Shell script
- [ ] Other (please specify): 


### Supporting software requirements

#### Version of primary software used
-->
<!--
(e.g., R version 3.6.0)
-->
<!--

#### Libraries and dependencies used by the code


Include version numbers (e.g., version numbers for any R or Python packages used)
-->

### Supporting system/hardware requirements (optional)

Simulations were conducted using parallel computing on 
the laboratory supercomputer, where 4 nodes and 36 cores are used. 


<!--
OPTIONAL: System/hardware requirements including operating system with version number, access to cluster, GPUs, etc.


### Parallelization used

- [ ] No parallel code used
- [ ] Multi-core parallelization on a single machine/node
    - Number of cores used: 
- [ ] Multi-machine/multi-node parallelization 
    - Number of nodes and cores used: 

### License

- [ ] MIT License (default)
- [ ] BSD 
- [ ] GPL v3.0
- [ ] Creative Commons
- [ ] Other: (please specify)


### Additional information (optional)

<!--
OPTIONAL: By default, submitted code will be published on the JASA GitHub repository (http://github.com/JASA-ACS) as well as in the supplementary material. Authors are encouraged to also make their code available in a public code repository, such as on GitHub, GitLab, or BitBucket. If relevant, please provide unique identifier/DOI/version information (e.g., a Git commit ID, branch, release, or tag). If the code and workflow are provided together, this section may be omitted, with information provided in the "Location" section below.
-->

# 3: Reproducibility 

<!--
The materials provided should provide a straightforward way for reviewers and readers to reproduce analyses with as few steps as possible. 
-->

All tables and figures in the paper can be reproduced using the $R$ code provided in the Github repository: https://github.com/Nago2020/Jasacc. There are two subfolders in this repository: **SimulationStudy** and **DataApplication**. In the subfolder **SimulationStudy**, all workflow information to reproduce tables and figures in the simulation study is contained in the **Master_reproducibility_simulations.R** script. Running this script on a desktop or laptop is very time-consuming, even after parallel processing is done via **foreach()**. 


The workflow information is provided in the **Real_data_analysis_results.R** script in the **DataApplication** subfolder for reproducing tables and figures in the data application section. **Real_data_analysis_results.R** requires approximately 5-6 hours to reproduce the results in Tables 4-5 and Figure 5 and the goodness-of-fit test results in Table 8.   


<!--
## Scope

The provided workflow reproduces:

- [ ] Any numbers provided in text in the paper
- [ ] The computational method(s) presented in the paper (i.e., code is provided that implements the method(s))
- [ ] All tables and figures in the paper
- [ ] Selected tables and figures in the paper, as explained and justified below:

-->

<!--
## Workflow

### Location

The workflow is available:
-->

<!--
Check all that apply, and in the case of a Git repository include unique identifier, such as specific commit ID, branch, release, or tag.
-->
<!--
- [ ] As part of the paper’s supplementary material.
- [ ] In this Git repository:
- [ ] Other (please specify):
-->
<!--
Indicate where the materials (generally including the code, unless in a separate location and indicated in the previous section) are available. We strongly encourage authors to place their materials (but not large datasets) in a Git repository hosted on a site such as GitHub, GitLab, or BitBucket. If the repository is private during the review process, please indicate the location where it will be available publicly upon publication, and also include the materials as a zip file (e.g., obtained directly from the Git hosting site) as supplementary materials.
-->

<!--
### Format(s)


Check all that appl
- [ ] Single master code file 
- [ ] Wrapper (shell) script(s)
- [ ] Self-contained R Markdown file, Jupyter notebook, or other literate programming approach
- [ ] Text file (e.g., a readme-style file) that documents workflow
- [ ] Makefile
- [ ] Other (more detail in *Instructions* below)

### Instructions

-->
<!--
Describe how to use the materials provided to reproduce analyses in the manuscript. Additional details can be provided in file(s) accompanying the reproducibility materials. If no workflow is provided, please state this and say why (e.g., if the paper contains no computational work).
-->

<!--
### Expected run-time

Approximate time needed to reproduce the analyses on a standard desktop machine:

- [ ] < 1 minute
- [ ] 1-10 minutes
- [ ] 10-60 minutes
- [ ] 1-8 hours
- [ ] > 8 hours
- [ ] Not feasible to run on a desktop machine, as described here:

### Additional information (optional)

<!--
OPTIONAL: Additional documentation provided (e.g., R package vignettes, demos or other examples) that show how to use the provided code/software in other settings.
-->
<!--
# Notes (optional)

<!--
OPTIONAL: Any other relevant information not covered on this form. If reproducibility materials are not publicly available at the time of submission, please provide information here on how the reviewers can view the materials.
-->
