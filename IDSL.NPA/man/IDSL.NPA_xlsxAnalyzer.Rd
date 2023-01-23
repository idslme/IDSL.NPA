\name{IDSL.NPA_xlsxAnalyzer}
\alias{IDSL.NPA_xlsxAnalyzer}
\title{
IDSL.NPA workflow xlsx Analyzer
}
\description{
This function processes the spreadsheet of the NPA parameters to ensure the parameter inputs are consistent with the requirements of the IDSL.NPA pipeline.
}
\usage{
IDSL.NPA_xlsxAnalyzer(spreadsheet)
}
\arguments{
  \item{spreadsheet}{
'Start' tab of the parameter spreadsheet
}
}
\value{
This function returns the NPA parameters to feed the `IDSL.NPA_workflow` function.
}
