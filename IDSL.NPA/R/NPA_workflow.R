NPA_workflow <- function(PARAM_NPA) {
  ##
  NPA0001 <- tolower(PARAM_NPA[which(PARAM_NPA[, 1] == 'NPA0001'), 2])
  NPA0002 <- tolower(PARAM_NPA[which(PARAM_NPA[, 1] == 'NPA0002'), 2])
  NPT <- as.numeric(PARAM_NPA[which(PARAM_NPA[, 1] == 'NPA0003'), 2])
  input_path_ms <- PARAM_NPA[which(PARAM_NPA[, 1] == 'NPA0004'), 2]
  ##
  output_address <- PARAM_NPA[which(PARAM_NPA[, 1] == 'NPA0008'), 2]
  FSA_dir.create(output_address, allowedUnlink = FALSE)
  opendir(output_address)
  ##
  ##############################################################################
  ## To create log record for IDSL.FSA
  initiation_time <- Sys.time()
  timeZone <- tryCatch(Sys.timezone(), warning = function(w) {"UTC"}, error = function(e) {"UTC"})
  .logFSA <- NULL
  .logFSA <<- paste0(output_address, "/logNPA_performance.txt")
  FSA_logRecorder(paste0(rep("", 100), collapse = "="))
  FSA_logRecorder(paste0("mzML/mzXML/netCDF:  ", input_path_ms))
  FSA_logRecorder(paste0("OUTPUT:  ", output_address))
  FSA_logRecorder(paste0(rep("", 100), collapse = "-"))
  FSA_logRecorder("Initiated Nominal Peak Analysis (NPA)!")
  FSA_logRecorder(paste0(as.character(initiation_time), " ", timeZone))
  FSA_logRecorder("", allowedPrinting = FALSE)
  FSA_logRecorder("", allowedPrinting = FALSE)
  FSA_logRecorder(paste0(PARAM_NPA[, 1], "\t", PARAM_NPA[, 2]),  allowedPrinting = FALSE)
  FSA_logRecorder(paste0(rep("", 100), collapse = "-"))
  ##
  ##############################################################################
  ##
  ref_xlsx_file <- as.character(PARAM_NPA[which(PARAM_NPA[, 1] == 'NPA0005'), 2])
  refMSPcreationCheck <- file.exists(ref_xlsx_file)
  ##
  ##############################################################################
  ##
  if (refMSPcreationCheck) {
    ##
    refNPAtable <- NPA_reference_xlsxAnalyzer(ref_xlsx_file, input_path_ms)[[1]]
    ##
    refMSindexList <- FSA_R.aggregate(refNPAtable$`Filename`)
    file_name_ms <- names(refMSindexList)
    mzRef <- as.numeric(refNPAtable$`MZ`)
    RTref <- as.numeric(refNPAtable$`RT`)
    refNPAtable$`MZ` <- NULL
    refNPAtable$`RT` <- NULL
    ##
    RTtoleranceRef <- as.numeric(PARAM_NPA[which(PARAM_NPA[, 1] == 'NPA0021'), 2])
    ##
    output_NPA_MSP <- paste0(output_address, "/NPA_REF_MSP")
    tryCatch(FSA_dir.create(output_NPA_MSP, allowedUnlink = TRUE), error = function(e) {message("ERROR!!! The .msp files inside the `NPA_REF_MSP` folder may influence the final output!")})
    ##
    str_xlsx_file <- strsplit(ref_xlsx_file, "/")[[1]]
    mspFileName <- str_xlsx_file[length(str_xlsx_file)]
    mspFileName <- gsub("[.]xlsx$|[.]csv$", ".msp", mspFileName, ignore.case = TRUE)
    mspFileName <- paste0("NPA_REF_MSP_", mspFileName)
    ##
  } else {
    refMSindexList <- NULL
    refNPAtable <- NULL
    massErrorRef <- 0
    RTtoleranceRef <- 0
    ##
    samples_string <- PARAM_NPA[which(PARAM_NPA[, 1] == 'NPA0006'), 2]
    if (tolower(samples_string) == "all") {
      file_name_ms <- dir(path = input_path_ms)
      file_name_ms <- file_name_ms[grep(pattern = ".mzML$|.mzXML$|.CDF$", file_name_ms, ignore.case = TRUE)]
    } else {
      file_name_ms <- strsplit(samples_string, ";")[[1]]
    }
    ##
    output_NPA_MSP <- paste0(output_address, "/NPA_MSP")
    FSA_dir.create(output_NPA_MSP, allowedUnlink = FALSE)
  }
  ##
  ##############################################################################
  ##
  IDSL.IPA::opendir(output_NPA_MSP)
  ##
  LMS <- length(file_name_ms)
  if (LMS == 0) {
    stop(FSA_logRecorder("EMPTY HRMS FOLDER!!!"))
  }
  ##
  ##############################################################################
  ##
  if (NPA0001 == "yes") {
    plotEICcheck <- if (tolower(PARAM_NPA[which(PARAM_NPA[, 1] == 'NPA0007'), 2]) == "yes") {TRUE} else {FALSE}
    if (plotEICcheck) {
      dev.offCheck <- TRUE
      while (dev.offCheck) {
        dev.offCheck <- tryCatch(dev.off(), error = function(e) {FALSE})
      }
      ##
      output_NPA_EICs_folder <- paste0(output_address, "/NPA_EIC")
      FSA_dir.create(output_NPA_EICs_folder, allowedUnlink = FALSE)
    } else {
      outputNPAeic <- NULL
      output_NPA_EICs_folder <- NULL
    }
    ##
    parallelizationMode <- tolower(PARAM_NPA[which(PARAM_NPA[, 1] == 'NPA0009'), 2])
    RTtolerance <- as.numeric(PARAM_NPA[which(PARAM_NPA[, 1] == 'NPA0010'), 2])
    minSNRbaseline <- as.numeric(PARAM_NPA[which(PARAM_NPA[, 1] == 'NPA0011'), 2])
    peakHeightThreshold <- as.numeric(PARAM_NPA[which(PARAM_NPA[, 1] == 'NPA0012'), 2])
    smoothingWindow <- as.numeric(PARAM_NPA[which(PARAM_NPA[, 1] == 'NPA0013'), 2])
    nSpline <- as.numeric(PARAM_NPA[which(PARAM_NPA[, 1] == 'NPA0014'), 2])
    topRatioPeakHeight <- as.numeric(PARAM_NPA[which(PARAM_NPA[, 1] == 'NPA0015'), 2])/100 # ratio of top peak percentage of chromatographic peaks to measure peak similarities (%)
    minIonRangeDifference <- as.numeric(PARAM_NPA[which(PARAM_NPA[, 1] == 'NPA0016'), 2])
    minNumNPApeaks <- as.numeric(PARAM_NPA[which(PARAM_NPA[, 1] == 'NPA0017'), 2])
    pearsonRHOthreshold <- as.numeric(PARAM_NPA[which(PARAM_NPA[, 1] == 'NPA0018'), 2])
    ##
    if (refMSPcreationCheck) {
      FSA_logRecorder("Individual `.msp` files are stored in the `NPA_REF_MSP` folder!")
    } else {
      FSA_logRecorder("Individual `.msp` files are stored in the `NPA_MSP` folder!")
    }
    ##
    ############################################################################
    ##
    NPA_workflow_call <- function(plotEICcheck, outputNPAeic, output_NPA_EICs_folder, input_path_ms, iMSfilename,
                                  smoothingWindow, peakHeightThreshold, minSNRbaseline, RTtolerance, nSpline,
                                  topRatioPeakHeight, minIonRangeDifference, minNumNPApeaks, pearsonRHOthreshold,
                                  NPT, refMSPcreationCheck, refMSindexList, refNPAtable, RTtoleranceRef, output_NPA_MSP) {
      ##
      if (plotEICcheck) {
        outputNPAeic <- paste0(output_NPA_EICs_folder, "/NPA_EICs_", iMSfilename)
      }
      ##
      NPA_peaklist <- NPA_fragmentationPeakDetection(input_path_ms, iMSfilename, smoothingWindow, peakHeightThreshold, minSNRbaseline,
                                                     RTtolerance, nSpline, topRatioPeakHeight, minIonRangeDifference, minNumNPApeaks,
                                                     pearsonRHOthreshold, outputNPAeic, number_processing_threads = NPT)
      ##
      if (NPA_peaklist[1, 1] != 0) {
        if (refMSPcreationCheck) {
          ## To update NPA peaklist to generate reference .msp files
          xRef <- refMSindexList[[iMSfilename]]
          listMatchedNPApeakIDs <- lapply(xRef, function(j) {
            xNPA <- which((NPA_peaklist[, 5] == mzRef[j]) & (abs(NPA_peaklist[, 3] - RTref[j]) <= RTtoleranceRef))
            ##
            if (length(xNPA) > 0) {
              xNPA1 <- which(NPA_peaklist[, 1] == NPA_peaklist[xNPA, 1])
              list(NPA_peaklist[xNPA1, 1], j)
            }
          })
          ##
          if (length(listMatchedNPApeakIDs) > 0) {
            xMatchedNPApeakIDs <- unique(do.call(c, lapply(listMatchedNPApeakIDs, function(j) {j[[1]][1]})))
            ##
            NPA_peaklist <- NPA_peaklist[(NPA_peaklist[, 1] %in% xMatchedNPApeakIDs), ]
            ##
            selectedIPApeaks_IDref <- do.call(rbind, lapply(listMatchedNPApeakIDs, function(j) {c(j[[1]][1], j[[2]])}))
            ##
            NPA_REF_MSP <- IDSL.NPA_referenceMSPgenerator(NPA_peaklist, refNPAtable, selectedIPApeaks_IDref)
            write.table(NPA_REF_MSP, file = paste0(output_NPA_MSP, "/NPA_REF_MSP_", iMSfilename, ".msp"), quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE)
            ##
          }
        } else {
          NPA_MSP <- IDSL.NPA_MSPgenerator(NPA_peaklist, number_processing_threads = NPT)
          write.table(NPA_MSP, file = paste0(output_NPA_MSP, "/NPA_MSP_", iMSfilename, ".msp"), quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE)
        }
      } else {
        FSA_logRecorder(paste0("No peak was detected for `", iMSfilename, "`!"))
      }
      ##
      return()
    }
    ##
    ##############################################################################
    ##
    if (NPT == 1 | parallelizationMode == "peakmode") {
      ##
      iCounter <- 0
      progressBARboundaries <- txtProgressBar(min = 0, max = LMS, initial = 0, style = 3)
      for (i in file_name_ms) {
        ##
        null_variable <- tryCatch(NPA_workflow_call(plotEICcheck, outputNPAeic, output_NPA_EICs_folder, input_path_ms, iMSfilename = i,
                                                    smoothingWindow, peakHeightThreshold, minSNRbaseline, RTtolerance, nSpline,
                                                    topRatioPeakHeight, minIonRangeDifference, minNumNPApeaks, pearsonRHOthreshold,
                                                    NPT, refMSPcreationCheck, refMSindexList, refNPAtable, RTtoleranceRef, output_NPA_MSP),
                                  error = function(e) {FSA_logRecorder(paste0("Problem with `", i,"`!"))})
        ##
        iCounter <- iCounter + 1
        setTxtProgressBar(progressBARboundaries, iCounter)
      }
      close(progressBARboundaries)
      ##
    } else if (parallelizationMode == "samplemode") {
      NPT0 <- NPT
      NPT <- 1
      ##
      osType <- Sys.info()[['sysname']]
      ##
      if (osType == "Linux") {
        ##
        null_variable <- mclapply(file_name_ms, function(i) {
          ##
          tryCatch(NPA_workflow_call(plotEICcheck, outputNPAeic, output_NPA_EICs_folder, input_path_ms, iMSfilename = i,
                                     smoothingWindow, peakHeightThreshold, minSNRbaseline, RTtolerance, nSpline,
                                     topRatioPeakHeight, minIonRangeDifference, minNumNPApeaks, pearsonRHOthreshold,
                                     NPT, refMSPcreationCheck, refMSindexList, refNPAtable, RTtoleranceRef, output_NPA_MSP),
                   error = function(e) {FSA_logRecorder(paste0("Problem with `", i,"`!"))})
        }, mc.cores = NPT0)
        ##
        closeAllConnections()
        ##
      } else if (osType == "Windows") {
        ##
        clust <- makeCluster(NPT0)
        registerDoParallel(clust)
        ##
        null_variable <- foreach(i = file_name_ms, .verbose = FALSE) %dopar% {
          ##
          tryCatch(NPA_workflow_call(plotEICcheck, outputNPAeic, output_NPA_EICs_folder, input_path_ms, iMSfilename = i,
                                     smoothingWindow, peakHeightThreshold, minSNRbaseline, RTtolerance, nSpline,
                                     topRatioPeakHeight, minIonRangeDifference, minNumNPApeaks, pearsonRHOthreshold,
                                     NPT, refMSPcreationCheck, refMSindexList, refNPAtable, RTtoleranceRef, output_NPA_MSP),
                   error = function(e) {FSA_logRecorder(paste0("Problem with `", i,"`!"))})
        }
        ##
        stopCluster(clust)
        ##
      }
      NPT <- NPT0
    }
    ##
    ############################################################################
    ##
    if (refMSPcreationCheck) {
      ref_msp_list <- dir(path = output_NPA_MSP, full.names = TRUE, pattern = ".msp$")
      if (length(ref_msp_list) > 0) {
        refMSP <- do.call(c, lapply(ref_msp_list, function(i) {
          readLines(i, warn = FALSE)
        }))
        ##
        write.table(refMSP, file = paste0(output_address, "/", mspFileName), quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE)
        FSA_logRecorder(paste0("The reference `", mspFileName, "` file was stored in the output directory!!!"))
      }
    }
  }
  ##
  ##############################################################################
  ##### Unique tag aggregation by spectra similarity across entire samples #####
  ##############################################################################
  ##
  if (NPA0002 == "yes") {
    plotSpectra <- if (tolower(PARAM_NPA[which(PARAM_NPA[, 1] == 'NPA0020'), 2]) == "yes") {TRUE} else {FALSE}
    allowedWeightedSpectralEntropy <- eval(parse(text = (PARAM_NPA[which(PARAM_NPA[, 1] == 'NPA0023'), 2])))
    minEntropySimilarity <- as.numeric(PARAM_NPA[which(PARAM_NPA[, 1] == 'NPA0024'), 2])
    ##
    if (refMSPcreationCheck) {
      if (file.exists(paste0(output_address, "/", mspFileName))) {
        FSA_logRecorder("Initiated detecting unique NPA variants!")
        ##
        aggregateBy <- PARAM_NPA[which(PARAM_NPA[, 1] == 'NPA0019'), 2]
        xAggregateBy <- which(colnames(refNPAtable) == aggregateBy)
        if (length(xAggregateBy) == 0) {
          aggregateBy <- "Name"
        }
        FSA_logRecorder(paste0("The meta-variable for aggregation is `", aggregateBy, "`!"))
        ##
        listSimilarMSPvariants <- FSA_uniqueMSPblockTagger(path = output_address, MSPfile = mspFileName, aggregateBy, massError = 0, RTtolerance = NA, minEntropySimilarity,
                                                           allowedNominalMass = TRUE, allowedWeightedSpectralEntropy, noiseRemovalRatio = 0, plotSpectra, number_processing_threads = NPT)
        FSA_logRecorder(paste0("Indices of similar MSP blocks for each compound are stored as `listSimilarMSPvariants.Rdata` in the `", output_address,"` folder!"))
        save(listSimilarMSPvariants, file = paste0(output_address, "/listSimilarMSPvariants.Rdata"))
        FSdb_address <- paste0(output_address, "/uniqueMSPtags_", gsub("[.]msp$|[.]Rdata$", ".Rdata", mspFileName, ignore.case = TRUE))
        FSA_logRecorder("Completed detecting unique NPA variants!")
      } else {
        FSdb_address <- ""
        FSA_logRecorder("No reference compound was detected!")
      }
      ##
    } else {
      ##
      RTtoleranceRef <- as.numeric(PARAM_NPA[which(PARAM_NPA[, 1] == 'NPA0021'), 2])
      minNPAdetectionFrequency <- floor(as.numeric(PARAM_NPA[which(PARAM_NPA[, 1] == 'NPA0022'), 2])*LMS/100)
      MSPfile_vector <- dir(path = output_NPA_MSP, pattern = ".msp$", ignore.case = TRUE)
      ##
      FSA_logRecorder("Initiated detecting unique NPA variants!")
      FSA_uniqueMSPblockTaggerUntargeted(path = output_NPA_MSP, MSPfile_vector, minNPAdetectionFrequency, minEntropySimilarity, massError = 0, massErrorPrecursor = NA, RTtoleranceRef,
                                         noiseRemovalRatio = 0, allowedNominalMass = TRUE, allowedWeightedSpectralEntropy, plotSpectra, number_processing_threads = NPT)
      FSdb_address <- paste0(output_NPA_MSP, "/UNIQUETAGS/uniqueMSPtagsUntargeted.Rdata")
      FSA_logRecorder("Completed detecting unique NPA variants!")
    }
  }
  ##
  ##############################################################################
  ##
  completion_time <- Sys.time()
  FSA_logRecorder(paste0(rep("", 100), collapse = "-"))
  required_time <- completion_time - initiation_time
  IPA_logRecorder(paste0("The required processing time was `", required_time, " ", attributes(required_time)$units, "`"))
  FSA_logRecorder(paste0(as.character(completion_time), " ", timeZone), allowedPrinting = FALSE)
  FSA_logRecorder("", allowedPrinting = FALSE)
  FSA_logRecorder("", allowedPrinting = FALSE)
  FSA_logRecorder("Completed the NPA analysis successfully!")
  FSA_logRecorder(paste0(rep("", 100), collapse = "="), allowedPrinting = FALSE)
  ##
  ##############################################################################
  ##
  return()
}
