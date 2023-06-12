IDSL.NPA_MSPgenerator <- function(NPA_peaklist, number_processing_threads = 1) {
  ##
  ID <- as.numeric(NPA_peaklist[, 1])
  x_diffID <- c(0, which(abs(diff(ID)) > 0), length(ID))
  ##
  call_MSPvector <- function(j) {
    ##
    x_ID <- seq((x_diffID[j] + 1), x_diffID[j + 1], 1)
    ID_j <-  NPA_peaklist[x_ID[1], 1]
    ##
    MSPid <- paste0("Name: NPApeakGrouping_ID_", ID_j, "_RT_", NPA_peaklist[x_ID[1], 3], "\n")
    ##
    MSPid <- paste0(MSPid, "NPA_ID: ", ID_j, "\n")
    ##
    MSPid <- paste0(MSPid, "Retention_Time: ", NPA_peaklist[x_ID[1], 3], "\n")
    ##
    MSPid <- paste0(MSPid, "SNRbaseline: ", paste0(NPA_peaklist[x_ID, 10], collapse = ","), "\n")
    ##
    MSPid <- paste0(MSPid, "Pearson_rho: ", paste0(NPA_peaklist[x_ID, 9], collapse = ","), "\n")
    ##
    MSPid <- paste0(MSPid, "MS_level: ", 1, "\n")
    ##
    IonMode <- if (as.numeric(NPA_peaklist[j, 8]) == 1) {"Positive"} else {"Negative"}
    ##
    MSPid <- paste0(MSPid, "Ion_mode: ", IonMode, "\n")
    ##
    MSPid <- paste0(MSPid, "Weighted_spectral_entropy_0noiseRemoval: ", NPA_peaklist[x_ID[1], 7], "\n")
    ##
    MSPid <- paste0(MSPid, "basePeakMZ: ", NPA_peaklist[x_ID[1], 5], "\n")
    MSPid <- paste0(MSPid, "basePeakIntensity: ", NPA_peaklist[x_ID[1], 6], "\n")
    ##
    MSPid <- paste0(MSPid, "Num Peaks: ", length(x_ID), "\n")
    ##
    MSPid_mz_int <- paste0(NPA_peaklist[x_ID, 5], " ", NPA_peaklist[x_ID, 6], "\n", collapse = "")
    ##
    paste0(MSPid, MSPid_mz_int, "\n")
  }
  ##############################################################################
  if (number_processing_threads == 1) {
    MSPvector <- do.call(c, lapply(1:(length(x_diffID) - 1), function(j) {
      call_MSPvector(j)
    }))
  } else {
    ##
    osType <- Sys.info()[['sysname']]
    ##
    if (osType == "Windows") {
      ##
      clust <- makeCluster(number_processing_threads)
      clusterExport(clust, setdiff(ls(), c("clust")), envir = environment())
      ##
      MSPvector <- do.call(c, parLapply(clust, 1:(length(x_diffID) - 1), function(i) {
        call_MSPvector(i)
      }))
      ##
      stopCluster(clust)
      ####
    } else {
      ##
      MSPvector <- do.call(c, mclapply(1:(length(x_diffID) - 1), function(j) {
        call_MSPvector(j)
      }, mc.cores = number_processing_threads))
      ##
      closeAllConnections()
      ##
    }
  }
  ##
  return(MSPvector)
}
