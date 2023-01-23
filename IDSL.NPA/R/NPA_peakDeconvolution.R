NPA_peakDeconvolution <- function(input_MS_path, MSfilename, smoothingWindow, peakHeightThreshold, minSNRbaseline, number_processing_threads = 1) {
  ############################ MS level = 1 ####################################
  p2l <- IDSL.MXP::peak2list(input_MS_path, MSfilename)
  scanTable <- p2l[["scanTable"]] # this gets table of details for each spectra
  spectraList <- p2l[["spectraList"]] # this gets the spectra values
  p2l <- NULL
  #
  x_MS <- which(scanTable$peaksCount > 0 & scanTable$msLevel == 1) # peaks from soft ionization channel ## some files may not have data in the column re-calibration period.
  spectraList <- spectraList[x_MS]
  scanTable <- scanTable[x_MS, ]
  ##
  RetentionTime <- scanTable$retentionTime
  nRT <- length(RetentionTime)       # n_RT is the maximum number of scan number
  ##############################################################################
  nRTSeq <- seq(1, nRT, 1)
  ##
  call_mzIntSCN <- function(j) {
    slj <- spectraList[[j]]
    cbind(slj, rep(j, dim(slj)[1]))
  }
  ##
  call_rawEIC <- function(j) {
    rawEICj <- rep(0, nRT)
    ##
    x_rt <- mzIntSCN[(index_xic[j] + 1):index_xic[j + 1], 3]
    rawEICj[x_rt] <- mzIntSCN[(index_xic[j] + 1):index_xic[j + 1], 2]
    ##
    rawEICj
  }
  ##
  call_smoothEIC <- function(j) {
    jChromatogram <- data.frame(cbind(nRTSeq, rawEIC[j, ]))
    colnames(jChromatogram) <- c("scan_number", "smooth_chrom")
    loess_SZC <- loess(smooth_chrom ~ scan_number, data = jChromatogram, span = smoothingWindow/nRT, control = loess.control(surface = "direct"))
    smoothJchromatogram <- predict(loess_SZC)
    x_neg <- which(smoothJchromatogram < 0)
    smoothJchromatogram[x_neg] <- 0
    ##
    smoothJchromatogram
  }
  ##
  call_peaklistNPA <- function(j) {
    x_apex <- which(matrixSegment[j, ] == +1)
    x_boundary <- which(matrixSegment[j, ] == -1)
    ##
    do.call(rbind, lapply(x_apex, function(k) {
      ##
      height <- rawEIC[j, k]
      if (height >= peakHeightThreshold) {
        snr <- height/baseline[j, k]
        if (!is.nan(snr)) {
          if (snr >= minSNRbaseline) {
            k11 <- which(x_boundary < k)
            k1 <- k11[length(k11)]
            k2 <- which(x_boundary > k)[1]
            ##
            c(j, k, mzEIC[j], x_boundary[k1], RetentionTime[k], x_boundary[k2], height, snr)
          }
        }
      }
    }))
  }
  ##############################################################################
  if (number_processing_threads == 1) {
    ##
    mzIntSCN <- do.call(rbind, lapply(1:nRT, function(j) {
      call_mzIntSCN(j)
    }))
    spectraList <- NULL
    ##
    mzIntSCN[, 1] <- round(mzIntSCN[, 1], 0)
    ##
    mzIntSCN <- mzIntSCN[order(mzIntSCN[, 1], decreasing = FALSE), ]
    ##
    index_xic <- c(0, which(diff(mzIntSCN[, 1]) > 0), dim(mzIntSCN)[1])
    Lindex_xic <- (length(index_xic) - 1)
    ##
    mzEIC <- do.call(c, lapply(1:Lindex_xic, function(j) {
      mzIntSCN[(index_xic[j] + 1), 1]
    }))
    ##
    rawEIC <- do.call(rbind, lapply(1:Lindex_xic, function(j) {
      call_rawEIC(j)
    }))
    mzIntSCN <- NULL
    ##
    smoothEIC <- do.call(rbind, lapply(1:Lindex_xic, function(j) {
      call_smoothEIC(j)
    }))
    ##
    matrixSegment <- do.call(rbind, lapply(1:Lindex_xic, function(j) {
      islocaloptimum(smoothEIC[j, ])
    }))
    ##
    baseline <- do.call(rbind, lapply(1:Lindex_xic, function(j) {
      Seg <- which(matrixSegment[j, ] == -1)
      IPA_baselineDeveloper(Seg, smoothEIC[j, ])
    }))
    ##
    peaklistNPA <- do.call(rbind, lapply(1:Lindex_xic, function(j) {
      call_peaklistNPA(j)
    }))
    ##
  } else {
    ##
    osType <- Sys.info()[['sysname']]
    ##
    if (osType == "Linux") {
      ##
      mzIntSCN <- do.call(rbind, mclapply(1:nRT, function(j) {
        call_mzIntSCN(j)
      }, mc.cores = number_processing_threads))
      spectraList <- NULL
      ##
      mzIntSCN[, 1] <- round(mzIntSCN[, 1], 0)
      ##
      mzIntSCN <- mzIntSCN[order(mzIntSCN[, 1], decreasing = FALSE), ]
      ##
      index_xic <- c(0, which(diff(mzIntSCN[, 1]) > 0), dim(mzIntSCN)[1])
      Lindex_xic <- (length(index_xic) - 1)
      ##
      mzEIC <- do.call(c, mclapply(1:Lindex_xic, function(j) {
        mzIntSCN[(index_xic[j] + 1), 1]
      }, mc.cores = number_processing_threads))
      ##
      rawEIC <- do.call(rbind, mclapply(1:Lindex_xic, function(j) {
        call_rawEIC(j)
      }, mc.cores = number_processing_threads))
      mzIntSCN <- NULL
      ##
      smoothEIC <- do.call(rbind, mclapply(1:Lindex_xic, function(j) {
        call_smoothEIC(j)
      }, mc.cores = number_processing_threads))
      ##
      matrixSegment <- do.call(rbind, mclapply(1:Lindex_xic, function(j) {
        islocaloptimum(smoothEIC[j, ])
      }, mc.cores = number_processing_threads))
      ##
      baseline <- do.call(rbind, mclapply(1:Lindex_xic, function(j) {
        Seg <- which(matrixSegment[j, ] == -1)
        IPA_baselineDeveloper(Seg, smoothEIC[j, ])
      }, mc.cores = number_processing_threads))
      ##
      peaklistNPA <- do.call(rbind, mclapply(1:Lindex_xic, function(j) {
        call_peaklistNPA(j)
      }, mc.cores = number_processing_threads))
      ##
      closeAllConnections()
      ##
    } else if (osType == "Windows") {
      ##
      clust <- makeCluster(number_processing_threads)
      registerDoParallel(clust)
      ##
      mzIntSCN <- foreach(j = 1:nRT, .combine = 'rbind', .verbose = FALSE) %dopar% {
        call_mzIntSCN(j)
      }
      spectraList <- NULL
      ##
      mzIntSCN[, 1] <- round(mzIntSCN[, 1], 0)
      ##
      mzIntSCN <- mzIntSCN[order(mzIntSCN[, 1], decreasing = FALSE), ]
      ##
      index_xic <- c(0, which(diff(mzIntSCN[, 1]) > 0), dim(mzIntSCN)[1])
      Lindex_xic <- (length(index_xic) - 1)
      ##
      mzEIC <- foreach(j = 1:Lindex_xic, .combine = 'c', .verbose = FALSE) %dopar% {
        mzIntSCN[(index_xic[j] + 1), 1]
      }
      ##
      rawEIC <- foreach(j = 1:Lindex_xic, .combine = 'rbind', .verbose = FALSE) %dopar% {
        call_rawEIC(j)
      }
      mzIntSCN <- NULL
      ##
      smoothEIC <- foreach(j = 1:Lindex_xic, .combine = 'rbind', .verbose = FALSE) %dopar% {
        call_smoothEIC(j)
      }
      ##
      matrixSegment <- foreach(j = 1:Lindex_xic, .combine = 'rbind', .verbose = FALSE) %dopar% {
        islocaloptimum(smoothEIC[j, ])
      }
      ##
      baseline <- foreach(j = 1:Lindex_xic, .combine = 'rbind', .verbose = FALSE) %dopar% {
        Seg <- which(matrixSegment[j, ] == -1)
        IPA_baselineDeveloper(Seg, smoothEIC[j, ])
      }
      ##
      peaklistNPA <- foreach(j = 1:Lindex_xic, .combine = 'rbind', .verbose = FALSE) %dopar% {
        call_peaklistNPA(j)
      }
      ##
      stopCluster(clust)
    }
  }
  ##
  outputList <- list(scanTable, rawEIC, smoothEIC, peaklistNPA)
  names(outputList) <- c("scanTable", "rawEIC", "smoothEIC", "peaklistNPA")
  ##
  return(outputList)
}
