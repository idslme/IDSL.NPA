\name{NPA_fragmentationPeakDetection}
\alias{NPA_fragmentationPeakDetection}
\title{
NPA Fragmentation Peaks Detection
}
\description{
This function detects fragmentation peaks for the NPA analysis
}
\usage{
NPA_fragmentationPeakDetection(input_MS_path, MSfilename, smoothingWindow,
peakHeightThreshold, minSNRbaseline, RTtolerance, nSpline, topRatioPeakHeight,
minIonRangeDifference, minNumNPApeaks, pearsonRHOthreshold, outputNPAeic = NULL,
number_processing_threads = 1)
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
  \item{RTtolerance}{
retention time tolerance to detect common peaks
}
  \item{nSpline}{
number of points for further smoothing using a cubic spline smoothing method to add more points to calculate Pearson correlation rho values
}
  \item{topRatioPeakHeight}{
The top percentage of the chromatographic peak to calculate Pearson correlation rho values
}
  \item{minIonRangeDifference}{
Minimum distance (Da) between lowest and highest m/z to prevent clustering isotopic envelopes
}
  \item{minNumNPApeaks}{
Minumum number of ions in a NPA cluster
}
  \item{pearsonRHOthreshold}{
Minimum threshold for Pearson correlation rho values
}
  \item{outputNPAeic}{
When `NULL` NPA EICs are not plotted. `outputNPAeic` represents an address to save NPA EICs figures.
}
  \item{number_processing_threads}{
Number of processing threads for multi-threaded processing
}
}
\value{
A dataframe peaklist of co-detected NPA analysis.
}
