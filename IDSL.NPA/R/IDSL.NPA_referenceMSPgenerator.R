IDSL.NPA_referenceMSPgenerator <- function(NPA_peaklist, refNPAtable, selectedPeaks_IDref) {
  ##
  x_1 <- which(NPA_peaklist$Ion_mode == "1")
  NPA_peaklist$Ion_mode[x_1] <- "Positive"
  ##
  x_0 <- which(NPA_peaklist$Ion_mode == "0")
  NPA_peaklist$Ion_mode[x_0] <- "Negative"
  ##
  ##############################################################################
  ##
  NPA_peaklist <- do.call(rbind, lapply(1:nrow(selectedPeaks_IDref), function(j) {
    x_j <- which(NPA_peaklist[, 1] == selectedPeaks_IDref[j, 1])
    l_x_j <- length(x_j)
    if (l_x_j > 0) {
      do.call(rbind, lapply(x_j, function(rowx_j) {
        cbind(NPA_peaklist[rowx_j, ], refNPAtable[selectedPeaks_IDref[j, 2], ])
      }))
    }
  }))
  ##
  ##############################################################################
  ##
  NPAcolNames <- colnames(NPA_peaklist)
  x_IDcol <- which(NPAcolNames == 'PeakID')
  x_mzCol <- which(NPAcolNames == 'mz_co_detected')
  x_intCol <- which(NPAcolNames == 'Int_co_detected')
  x_nameCol <- which(NPAcolNames == 'Name')
  ##
  if ((length(x_IDcol) != 1) | (length(x_mzCol) != 1) | (length(x_intCol) != 1) | (length(x_nameCol) != 1)) {
    stop("The NPA_peaklist file should have only one column for the following headers (case-senstive):
    'PeakID'
    'mz_co_detected'
    'Int_co_detected'
    'Name'")
  }
  ##
  NPAcolumns <- setdiff(NPAcolNames, c('PeakID', 'mz_co_detected', 'Int_co_detected', 'Name'))
  if (length(NPAcolumns) > 0) {
    NPAcolumnsCheck <- TRUE
    x_NPAcolumns <- do.call(c, lapply(NPAcolumns, function(i) {
      which(NPAcolNames == i)
    }))
  } else {
    NPAcolumnsCheck <- FALSE
  }
  ##
  PeakID <- as.numeric(NPA_peaklist[, x_IDcol])
  x_diffID <- c(0, which(abs(diff(PeakID)) > 0), length(PeakID))
  ##
  MSPvector <- do.call(c, lapply(1:(length(x_diffID) - 1), function(i) {
    x_ID <- seq((x_diffID[i] + 1), x_diffID[i + 1], 1)
    ID_i <-  x_diffID[i] + 1
    ##
    MSP1 <- paste0("Name: ", NPA_peaklist[ID_i, x_nameCol], "\n")
    MSP1 <- paste0(MSP1, "MSP_mode: NPA\n")
    MSP1 <- paste0(MSP1, "MS_level: 1\n")
    ##
    MSP1 <- paste0(MSP1, "basePeakMZ: ", NPA_peaklist[x_ID[1], x_mzCol], "\n")
    MSP1 <- paste0(MSP1, "basePeakIntensity: ", NPA_peaklist[x_ID[1], x_intCol], "\n")
    ##
    if (NPAcolumnsCheck) {
      MSPid <- do.call(paste0, lapply(x_NPAcolumns, function(j) {
        paste0(NPAcolNames[j], ": ", NPA_peaklist[ID_i, j], "\n")
      }))
    } else {
      MSPid <- NULL
    }
    ##
    MSPid <- paste0(MSPid, "Num Peaks: ", length(x_ID), "\n")
    ##
    MSPid_mz_int <- paste0(NPA_peaklist[x_ID, x_mzCol], " ", NPA_peaklist[x_ID, x_intCol], "\n", collapse = "")
    ##
    paste0(MSP1, MSPid, MSPid_mz_int, "\n")
  }))
  ##
  return(MSPvector)
}
