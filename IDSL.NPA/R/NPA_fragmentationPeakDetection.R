NPA_fragmentationPeakDetection <- function(input_MS_path, MSfilename, smoothingWindow, peakHeightThreshold, minSNRbaseline,
                                           RTtolerance, nSpline, topRatioPeakHeight, minIonRangeDifference, minNumNPApeaks,
                                           pearsonRHOthreshold, outputNPAeic = NULL, number_processing_threads = 1) {
  ##
  minNumNPApeaks <- minNumNPApeaks - 1
  ##
  ##############################################################################
  if (is.null(outputNPAeic)) {
    plotEICcheck <- FALSE
  } else {
    ##
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    ##
    plotEICcheck <- TRUE
    ##
    FSA_dir.create(outputNPAeic, allowedUnlink = FALSE)
  }
  ############################ MS level = 1 ####################################
  deconvolutionList <- NPA_peakDeconvolution(input_MS_path, MSfilename, smoothingWindow, peakHeightThreshold, minSNRbaseline, number_processing_threads)
  scanTable <- deconvolutionList[["scanTable"]]
  rawEIC <- deconvolutionList[["rawEIC"]]
  smoothEIC <- deconvolutionList[["smoothEIC"]]
  peaklistNPA <- deconvolutionList[["peaklistNPA"]]
  deconvolutionList <- NULL
  ##
  RetentionTime <- scanTable$retentionTime
  MS1_polarity <- as.matrix(scanTable$polarity)
  scanTable <- NULL
  ##############################################################################
  counterNPAblock <- 0
  if (!is.null(peaklistNPA)) {
    nPeaks <- dim(peaklistNPA)[1]
    peaklistNPA <- peaklistNPA[order(peaklistNPA[, 8], decreasing = TRUE), ]
    idNPA <- peaklistNPA[, 1]
    apexID <- peaklistNPA[, 2]
    leftBoundary <- peaklistNPA[, 4]
    rightBoundary <- peaklistNPA[, 6]
    mzNPA <- peaklistNPA[, 3]
    rtNPA <- peaklistNPA[, 5]
    intNPA <- peaklistNPA[, 7]
    SNRbaseline <- peaklistNPA[, 8]
    ############################################################################
    j <- 1
    listNPApeaklist <- vector(mode = "list", nPeaks)
    while (j < nPeaks) {
      if (idNPA[j] != 0) {
        ########################################################################
        x_fragment <- which((abs(rtNPA - rtNPA[j]) <= RTtolerance) & (idNPA != 0))
        iMzNPAions <- mzNPA[x_fragment]
        ionRangeDifference <- max(iMzNPAions) - min(iMzNPAions)
        if (ionRangeDifference >= minIonRangeDifference) {
          ######################################################################
          x_fragment <- setdiff(x_fragment, j)
          ##
          L_fragment <- length(x_fragment)
          if (L_fragment >= minNumNPApeaks) {
            ##
            RT_chrom_precursor <- RetentionTime[leftBoundary[j]:rightBoundary[j]]
            Int_chrom_precursor <- smoothEIC[idNPA[j], leftBoundary[j]:rightBoundary[j]]
            ##
            W_precursor <- spline(RT_chrom_precursor, Int_chrom_precursor , n = nSpline, method = "fmm",
                                  xmin = RT_chrom_precursor[1], xmax = RT_chrom_precursor[length(RT_chrom_precursor)], ties = mean) # To smooth the curve for derivative calculations
            RT_spline_precursor <- W_precursor[[1]]
            Int_spline_precursor <- W_precursor[[2]]
            #
            x_topRatioPeakHeight <- which(Int_spline_precursor/max(Int_spline_precursor) >= (1 - topRatioPeakHeight))
            RT_spline_precursor <- RT_spline_precursor[x_topRatioPeakHeight]
            Int_spline_precursor <- Int_spline_precursor[x_topRatioPeakHeight]
            ####################################################################
            NPA_EICs <- lapply(x_fragment, function(k) {
              ##
              RT_chrom_fragment <- RetentionTime[leftBoundary[k]:rightBoundary[k]]
              Int_chrom_fragment <- smoothEIC[idNPA[k], leftBoundary[k]:rightBoundary[k]]
              ##
              height_fragment <- rawEIC[idNPA[k], apexID[k]]
              if (height_fragment > 0) {
                ##
                Int_spline_fragment <- approx(RT_chrom_fragment, Int_chrom_fragment, RT_spline_precursor, method = "linear", 0, 0, rule = 2, f = 0, ties = mean)[[2]]
                #
                pearsonRHO <- suppressWarnings(cor(Int_spline_precursor, Int_spline_fragment, method = "pearson"))
                ##
                if (!is.na(pearsonRHO)) {
                  if (pearsonRHO >= pearsonRHOthreshold) {
                    if (plotEICcheck) {
                      L_chrom_fragment <- length(Int_chrom_fragment)
                      NPAEICdata <- cbind(rep(k, L_chrom_fragment), RT_chrom_fragment, Int_chrom_fragment)
                    } else {
                      NPAEICdata <- NULL
                    }
                    ##
                    NPA_fragments <- c(mzNPA[k], height_fragment, pearsonRHO, SNRbaseline[k])
                    ##
                    list(NPAEICdata, NPA_fragments, k)
                  }
                }
              }
            })
            ####################################################################
            x_fragment_npa <- do.call(c, lapply(1:L_fragment, function(k) {
              if (!is.null(NPA_EICs[[k]])) {
                k
              }
            }))
            ##
            if (length(x_fragment_npa) >= minNumNPApeaks) {
              NPA_fragments <- do.call(rbind, lapply(x_fragment_npa, function(k) {
                NPA_EICs[[k]][[2]]
              }))
              ##################################################################
              peakIDjk <- c(j, do.call(c, lapply(x_fragment_npa, function(k) {
                NPA_EICs[[k]][[3]]
              })))
              ##################################################################
              NPAmz <- c(mzNPA[j], intNPA[j], 1, SNRbaseline[j])
              NPAfragmentationList <- rbind(NPAmz, NPA_fragments)
              ##
              ionRangeDifference <- max(NPAfragmentationList[, 1]) - min(NPAfragmentationList[, 1])
              if (ionRangeDifference >= minIonRangeDifference) {
                ##
                NPAfragmentationList <- NPAfragmentationList[order(NPAfragmentationList[, 2], decreasing = TRUE), ]
                ##
                spectralEntropy <- round(spectral_entropy_calculator(NPAfragmentationList[, 1:2], allowedWeightedSpectralEntropy = TRUE, noiseRemovalRatio = 1e-16)[[1]], 5)
                ##
                counterNPAblock <- counterNPAblock + 1
                ##
                L_npa_fragments <- nrow(NPAfragmentationList)
                NPA_list <- cbind(rep(counterNPAblock, L_npa_fragments), rep(mzNPA[j], L_npa_fragments), rep(rtNPA[j], L_npa_fragments), rep(intNPA[j], L_npa_fragments),
                                  NPAfragmentationList[, 1:2], rep(spectralEntropy, L_npa_fragments), rep(MS1_polarity[apexID[j]], L_npa_fragments), NPAfragmentationList[, 3:4])
                ##
                idNPA[peakIDjk] <- 0
                ##
                listNPApeaklist[[counterNPAblock]] <- NPA_list
                ################################################################
                if (plotEICcheck) {
                  NPAEICdata <- do.call(rbind, lapply(x_fragment_npa, function(k) {
                    NPA_EICs[[k]][[1]]
                  }))
                  ##
                  yMaxLimPlot <- max(c(Int_chrom_precursor, NPAEICdata[, 3]))
                  ##
                  xLinesDiff <- c(0, which(abs(diff(NPAEICdata[, 1])) > 0), nrow(NPAEICdata))
                  nLines <- length(xLinesDiff)
                  nLines1 <- nLines - 1
                  ##
                  legText <- rep("", nLines)
                  legText[1] <- paste0("* m/z = ", mzNPA[j])
                  colors <- c("black", rainbow(nLines1, alpha = 1))
                  legCol <- rep("", nLines)
                  legCol[1] <- colors[1]
                  ##
                  orderMSP <- order(NPA_fragments[, 2], decreasing = TRUE)
                  orderLines <- do.call(rbind, lapply(2:nLines, function(p) {
                    c((xLinesDiff[p - 1] + 1), xLinesDiff[p])
                  }))
                  orderLines <- matrix(orderLines[orderMSP, ], ncol = 2)
                  ##
                  alignedEICfilename <- paste0(outputNPAeic, "/NPApeakGrouping_ID_", counterNPAblock, "_MZ_", mzNPA[j], "_RT_", rtNPA[j], ".png")
                  png(alignedEICfilename, width = 16, height = 8, units = "in", res = 100)
                  ##
                  par(mar = c(5.1, 4.1, 4.1, 13.8))
                  plot(RT_chrom_precursor, Int_chrom_precursor, type = "l", ylim = c(0, yMaxLimPlot*1.01), lwd = 4, col = colors[1], cex = 4, xlab = "", ylab = "")
                  ##
                  pCounter <- 1
                  for (p in 1:nLines1) {
                    pCounter <- pCounter + 1
                    ##
                    xLines <- seq(orderLines[p, 1], orderLines[p, 2], 1)
                    ##
                    lines(NPAEICdata[xLines, 2], NPAEICdata[xLines, 3], lwd = 2, col = colors[pCounter], cex = 4)
                    ##
                    legText[pCounter] <- paste0("m/z = ", mzNPA[NPAEICdata[xLines[1], 1]])
                    legCol[pCounter] <- colors[pCounter]
                  }
                  ##
                  mtext(text = paste0("S = ", spectralEntropy), side = 3, adj = 1, line = 0.25, cex = 1.0)
                  mtext("Retention time (min)", side = 1, adj = 0.5, line = 2, cex = 1.35)
                  mtext("Intensity", side = 2, adj = 0.5, line = 2, cex = 1.35)
                  mtext(MSfilename, side = 3, adj = 0, line = 0.25, cex = 1.4)
                  legend(x = "topright", inset = c(-0.22, 0), legend = legText, lwd = c(4, rep(2, nLines1)), cex = 1.125, bty = "n",
                         col = legCol, seg.len = 1, x.intersp = 0.5, y.intersp = 0.9, xpd = TRUE)
                  ##
                  dev.off()
                  ##############################################################
                }
              }
            }
          }
        }
      }
      ##
      x_j <- which(idNPA[(j + 1):nPeaks] != 0)
      if (length(x_j) > 0) {
        j <- j + x_j[1]
      } else {
        j <- nPeaks
      }
    }
  }
  ##############################################################################
  if (counterNPAblock > 0) {
    ##
    NPA_peaklist <- do.call(rbind, lapply(listNPApeaklist, function(j) {j}))
    ##
    NPA_peaklist <- data.frame(matrix(NPA_peaklist, ncol = 10))
    NPA_peaklist[, 1] <- as.numeric(NPA_peaklist[, 1])
    NPA_peaklist[, 5] <- round(as.numeric(NPA_peaklist[, 5]), 5)
    NPA_peaklist[, 6] <- round(as.numeric(NPA_peaklist[, 6]), 0)
    NPA_peaklist[, 9] <- round(as.numeric(NPA_peaklist[, 9]), 2)
    NPA_peaklist[, 10] <- round(as.numeric(NPA_peaklist[, 10]), 0)
    ##
  } else {
    NPA_peaklist <- data.frame(matrix(rep(0, 10), nrow = 1))
  }
  rownames(NPA_peaklist) <- NULL
  colnames(NPA_peaklist) <- c("PeakID", "mz12CIPA", "RTIPA", "IntIPA", "mz_co_detected", "Int_co_detected", "spectral_entropy", "Ion_mode", "rho", "SNRbaseline")
  ##
  return(NPA_peaklist)
}
