\name{NPA_peakDeconvolution}
\alias{NPA_peakDeconvolution}
\title{
NPA Peak Deconvolution
}
\description{
This function detects fragmentation peaks for the NPA analysis.
}
\usage{
NPA_peakDeconvolution(input_MS_path, MSfilename, smoothingWindow,
peakHeightThreshold, minSNRbaseline, number_processing_threads = 1)
}
\arguments{
  \item{input_MS_path}{
path to the MS files
}
  \item{MSfilename}{
MS file
}
  \item{smoothingWindow}{
number of scans for peak smoothing.
}
  \item{peakHeightThreshold}{
A minimum peak height threshold
}
  \item{minSNRbaseline}{
A minimum baseline S/N threshold
}
  \item{number_processing_threads}{
Number of processing threads for multi-threaded processing
}
}
\value{
A dataframe peaklist of co-detected DIA analysis.
}
