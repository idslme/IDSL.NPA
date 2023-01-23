\name{NPA_workflow}
\alias{NPA_workflow}
\title{
FSA NPA Workflow
}
\description{
This function runs the NPA analysis.
}
\usage{
NPA_workflow(PARAM_NPA)
}
\arguments{
  \item{PARAM_NPA}{
PARAM_NPA
}
}
\value{
This function stores .Rdata and .csv data from dataframe peaklist of co-detected NPA analysis.
}
\examples{
\donttest{
s_path <- system.file("extdata", package = "IDSL.NPA")
SSh1 <- paste0(s_path,"/NPA_parameters.xlsx")
## To see the results, use a known folder instead of the `tempdir()` command
temp_wd <- tempdir()
temp_wd_zip <- paste0(temp_wd, "/idsl.npa_testfile.zip")
spreadsheet <- readxl::read_xlsx(SSh1, sheet = "NPA")
PARAM_NPA <- cbind(spreadsheet[, 2], spreadsheet[, 4])
download.file(paste0("https://github.com/idslme/IDSL.NPA/blob/main/",
                     "NPA_educational_files/idsl.npa_testfile.zip?raw=true"),
              destfile = temp_wd_zip, mode = "wb")
unzip(temp_wd_zip, exdir = temp_wd)
PARAM_NPA[2, 2] <- "NO"
PARAM_NPA[4, 2] <- temp_wd
PARAM_NPA[8, 2] <- temp_wd
## To ensure `PARAM_NPA` is consistent with the `NPA_workflow`
PARAM_NPA <- NPA_xlsxAnalyzer(PARAM_NPA)
##
NPA_workflow(PARAM_NPA)
}
}
