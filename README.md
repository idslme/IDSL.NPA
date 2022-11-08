# IDSL.NPA<img src='NPA_educational_files/Figures/IDSL.NPA-logo.png' width="250px" align="right" />

<!-- badges: start -->
[![Maintainer](https://img.shields.io/badge/maintainer-Sadjad_Fakouri_Baygi-blue)](https://github.com/sajfb)
[![CRAN status](https://www.r-pkg.org/badges/version/IDSL.NPA)](https://cran.r-project.org/package=IDSL.NPA)
![](http://cranlogs.r-pkg.org/badges/IDSL.NPA?color=orange)
![](http://cranlogs.r-pkg.org/badges/grand-total/IDSL.NPA?color=brightgreen)
[![Dependencies](https://tinyverse.netlify.com/badge/IDSL.NPA)](https://cran.r-project.org/package=IDSL.NPA)
<!-- badges: end -->

**Nominal Peak Analysis (NPA)** by the [**Integrated Data Science Laboratory for Metabolomics and Exposomics (IDSL.ME)**](https://www.idsl.me/) is an R package to process nominal mass spectrometry data to create and annotate .msp files for untargeted MS/MS workflows.

	install.packages("IDSL.NPA") # IDSL.NPA package is set to release on CRAN by the end of 2022
	
## Workflow
You should download the [IDSL.NPA parameter spreadsheet](https://raw.githubusercontent.com/idslme/IDSL.NPA/main/NPA_parameters.xlsx) and select the parameters accordingly and then use this spreadsheet as the input for the IDSL.NPA workflow:	

	library(IDSL.NPA)
	IDSL.NPA_workflow("Address of the NPA parameter spreadsheet")

## Citation

pending...