\name{IDSL.NPA_MSPgenerator}
\alias{IDSL.NPA_MSPgenerator}
\title{
IDSL.NPA MSP Generator
}
\description{
This function creates standard .msp files that can also be used for Pepsearch.
}
\usage{
IDSL.NPA_MSPgenerator(NPA_peaklist, number_processing_threads = 1)
}
\arguments{
  \item{NPA_peaklist}{
A dataframe peaklist of co-detected peaks
}
  \item{number_processing_threads}{
Number of processing threads for multi-threaded processing
}
}
\value{
A string of standard .msp file
}
