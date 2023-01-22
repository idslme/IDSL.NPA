# IDSL.NPA<img src='NPA_educational_files/Figures/IDSL.NPA-logo.png' width="250px" align="right" />

<!-- badges: start -->
[![Maintainer](https://img.shields.io/badge/maintainer-Sadjad_Fakouri_Baygi-blue)](https://github.com/sajfb)
[![CRAN status](https://www.r-pkg.org/badges/version/IDSL.NPA)](https://cran.r-project.org/package=IDSL.NPA)
![](http://cranlogs.r-pkg.org/badges/IDSL.NPA?color=orange)
![](http://cranlogs.r-pkg.org/badges/grand-total/IDSL.NPA?color=brightgreen)
[![Dependencies](https://tinyverse.netlify.com/badge/IDSL.NPA)](https://cran.r-project.org/package=IDSL.NPA)
<!-- badges: end -->

**Nominal Peak Analysis (NPA)** by the [**Integrated Data Science Laboratory for Metabolomics and Exposomics (IDSL.ME)**](https://www.idsl.me/) is an R package to process nominal mass spectrometry data to create and annotate *.msp* files for untargeted MS/MS workflows.

## Table of Contents

- [Features of IDSL.NPA](https://github.com/idslme/IDSL.NPA#features-of-idslnpa)
- [Installation](https://github.com/idslme/IDSL.NPA#installation)
- [Workflow](https://github.com/idslme/IDSL.NPA#workflow)
- [Batch Example](https://github.com/idslme/IDSL.NPA#batch-example)
- [Citation](https://github.com/idslme/IDSL.NPA#citation)

## Installation

	install.packages("IDSL.NPA")
	
**Note:** In case you want to process **netCDF/CDF** mass spectrometry data, you should also install the [**RnetCDF**](https://CRAN.R-project.org/package=RNetCDF) package separately using the below command.

	install.packages("RNetCDF")
	
## Workflow
You should download the [IDSL.NPA parameter spreadsheet](https://raw.githubusercontent.com/idslme/IDSL.NPA/main/NPA_parameters.xlsx) and select the parameters accordingly and then use this spreadsheet as the input for the IDSL.NPA workflow:	

	library(IDSL.NPA)
	IDSL.NPA_workflow("Address of the NPA parameter spreadsheet")

## Batch Example

Follow these steps for a case study (n=300) [ST001154](https://www.metabolomicsworkbench.org/data/DRCCStudySummary.php?Mode=SetupRawDataDownload&StudyID=ST001154) which has Agilent 6890N data collected in EI positive mode.

1. Transfer mass spectrometry data with ***.cdf*** extension in a separate folder

2. The **Nominal Mass Analysis** requires 24 parameters distributed into 4 separate sections for a full scale analysis. For this study, use default parameter values presented in the [NPA parameter spreadsheet](https://raw.githubusercontent.com/idslme/IDSL.NPA/main/NPA_parameters.xlsx). Next, provide information for 
	
	2.1. Select **YES** for **PARAM0001** in the `Start` tab to only process **NPA** workflow.
	
	2.2. **NPA0004** for *Input data location (MS data)*
	
	2.3. **NPA0008** for *Output location (.msp files and EICs)*
	
	2.4. You may also increase the number of processing threads using **NPA0003** according to your computational power

3. Run this command in R/Rstudio console or terminal:

```
library(IDSL.NPA)
IDSL.NPA_workflow("Address of the NPA parameter spreadsheet")
```

4. You may parse the results at the address you provided for **NPA0008**.

## Citation

pending...