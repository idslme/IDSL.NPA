NPA_xlsxAnalyzer <- function(spreadsheet) {
  ##
  checkpoint_parameter <- FALSE
  ##
  if (typeof(spreadsheet) == "list") {
    if (ncol(spreadsheet) >= 4) {
      PARAM_NPA <- cbind(spreadsheet[, 2], spreadsheet[, 4])
      checkpoint_parameter <- TRUE
      ##
    } else if (ncol(spreadsheet) == 2) {
      PARAM_NPA <- spreadsheet
      checkpoint_parameter <- TRUE
      ##
    } else {
      FSA_message("The `NPA` spreadsheet tab was not produced properly!")
    }
  } else if (typeof(spreadsheet) == "character") {
    if (length(spreadsheet) == 1) {
      if (file.exists(spreadsheet)) {
        PARAM_NPA <- readxl::read_xlsx(spreadsheet, sheet = "NPA")
        PARAM_NPA <- cbind(PARAM_NPA[, 2], PARAM_NPA[, 4])
        checkpoint_parameter <- TRUE
      } else {
        FSA_message("The `NPA` spreadsheet tab not found! It should be an Excel file with .xlsx extention!")
      }
    } else {
      FSA_message("The `NPA` spreadsheet tab was not produced properly!")
    }
  } else {
    FSA_message("The `NPA` spreadsheet tab was not produced properly!")
  }
  ##############################################################################
  if (checkpoint_parameter) {
    ##
    x0001 <- which(PARAM_NPA[, 1] == 'NPA0001')
    if (length(x0001) == 0) {
      FSA_message("ERROR!!! Problem with NPA0001!")
      checkpoint_parameter <- FALSE
    } else {
      NPA0001 <- tolower(PARAM_NPA[x0001, 2])
      if (NPA0001 == "yes" | NPA0001 == "no") {
        PARAM_NPA[x0001, 2] <- NPA0001
      } else {
        FSA_message("ERROR!!! Problem with NPA0001!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    x0002 <- which(PARAM_NPA[, 1] == 'NPA0002')
    if (length(x0002) == 0) {
      FSA_message("ERROR!!! Problem with NPA0002!")
      checkpoint_parameter <- FALSE
    } else {
      NPA0002 <- tolower(PARAM_NPA[x0002, 2])
      if (NPA0002 == "yes" | NPA0002 == "no") {
        PARAM_NPA[x0002, 2] <- NPA0002
      } else {
        FSA_message("ERROR!!! Problem with NPA0002!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    number_processing_threads <- as.numeric(PARAM_NPA[which(PARAM_NPA[, 1] == 'NPA0003'), 2])
    if (length(number_processing_threads) == 0) {
      FSA_message("ERROR!!! Problem with NPA0003! This parameter should be a positive integer!")
      checkpoint_parameter <- FALSE
    } else {
      if (number_processing_threads >= 1) {
        if ((number_processing_threads %% 1) != 0) {
          FSA_message("ERROR!!! Problem with NPA0003! This parameter should be a positive integer!")
          checkpoint_parameter <- FALSE
        }
      } else {
        FSA_message("ERROR!!! Problem with NPA0003! This parameter should be at least 1 !")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    x0004 <- which(PARAM_NPA[, 1] == 'NPA0004')
    if (length(x0004) == 0) {
      FSA_message("ERROR!!! Problem with NPA0004!")
      checkpoint_parameter <- FALSE
    } else {
      input_path_ms <- PARAM_NPA[x0004, 2]
      input_path_ms <- gsub("\\", "/", input_path_ms, fixed = TRUE)
      PARAM_NPA[x0004, 2] <- input_path_ms
      if (!dir.exists(input_path_ms)) {
        FSA_message("ERROR!!! Problem with NPA0004! Please make sure the full path is provided!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    refMSPcreationCheck <- FALSE
    ref_xlsx_file <- as.character(PARAM_NPA[which(PARAM_NPA[, 1] == 'NPA0005'), 2])
    if (tolower(ref_xlsx_file) != "na") {
      refMSPcreationCheck <- TRUE
      ##
      listRefXlsxAnalyzer <- NPA_reference_xlsxAnalyzer(ref_xlsx_file, input_path_ms, PARAM = PARAM_NPA, PARAM_ID = 'NPA0005', checkpoint_parameter)
      ref_table <- listRefXlsxAnalyzer[[1]]
      PARAM_NPA <- listRefXlsxAnalyzer[[2]]
      checkpoint_parameter <- listRefXlsxAnalyzer[[3]]
      listRefXlsxAnalyzer <- NULL
      samples_string <- paste0(unique(ref_table$Filename), collapse = ";")
      ref_table <- NULL
    }
    ##
    LMS <- 0
    x0006 <- which(PARAM_NPA[, 1] == 'NPA0006')
    if (is.na(PARAM_NPA[x0006, 2])) {
      FSA_message("ERROR!!! Problem with NPA0006!")
      checkpoint_parameter <- FALSE
    } else {
      ##
      if (refMSPcreationCheck) {
        PARAM_NPA[x0006, 2] <- samples_string
      }
      ##
      if (tolower(PARAM_NPA[x0006, 2]) == "all") {
        file_name_ms <- dir(path = input_path_ms)
        file_name_ms <- file_name_ms[grep(pattern = ".mzML$|.mzXML$|.CDF$", file_name_ms, ignore.case = TRUE)]
        LMS <- length(file_name_ms)
        if (LMS == 0) {
          FSA_message("ERROR!!! Problem with NPA0006! No mzML/mzXML/CDF file was detected in the folder!")
        }
      } else {
        samples_string <- PARAM_NPA[x0006, 2]
        file_name_ms <- strsplit(samples_string, ";")[[1]]
        LMS <- length(file_name_ms)
        ndMS <- do.call(c, lapply(file_name_ms, function(i) {
          if (!file.exists(paste0(input_path_ms, "/", i))) {
            i
          }
        }))
        ##
        if (!is.null(ndMS)) {
          if (refMSPcreationCheck) {
            FSA_message("ERROR!!! The following file(s) can not be detected in the reference xlsx file (NPA0005) (case sensitive even for file extensions):")
          } else {
            FSA_message("ERROR!!! Problem with NPA0006! not detected the following file(s) (case sensitive even for file extensions):")
          }
          ##
          LMS <- LMS - length(ndMS)
          for (i in ndMS) {
            FSA_message(i)
          }
          checkpoint_parameter <- FALSE
        }
      }
    }
    ##
    x0007 <- which(PARAM_NPA[, 1] == 'NPA0007')
    if (length(x0007) == 0) {
      FSA_message("ERROR!!! Problem with NPA0007!")
      checkpoint_parameter <- FALSE
    } else {
      NPA0007 <- tolower(PARAM_NPA[x0007, 2])
      if (NPA0007 == "yes" | NPA0007 == "no") {
        PARAM_NPA[x0007, 2] <- NPA0007
      } else {
        FSA_message("ERROR!!! Problem with NPA0007!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    x0008 <- which(PARAM_NPA[, 1] == 'NPA0008')
    if (length(x0008) == 0) {
      FSA_message("ERROR!!! Problem with NPA0008!")
      checkpoint_parameter <- FALSE
    } else {
      output_path <- gsub("\\", "/", PARAM_NPA[x0008, 2], fixed = TRUE)
      PARAM_NPA[x0008, 2] <- output_path
      if (!dir.exists(output_path)) {
        tryCatch(dir.create(output_path, recursive = TRUE), warning = function(w){warning("Problem with NPA0008! R cannot create the folder!")})
        if (!dir.exists(output_path)) {
          checkpoint_parameter <- FALSE
        }
      }
    }
    ############################################################################
    ###### Chromatographic peak matching criteria to generate .msp files #######
    ############################################################################
    if (number_processing_threads > 1) {
      x0009 <- which(PARAM_NPA[, 1] == 'NPA0009')
      parallelizationMode <- PARAM_NPA[x0009, 2]
      if (is.na(parallelizationMode)) {
        FSA_message("ERROR!!! Problem with NPA0009!")
        checkpoint_parameter <- FALSE
      } else {
        parallelizationMode <- gsub(" ", "", tolower(parallelizationMode))
        if (parallelizationMode == "samplemode" | parallelizationMode == "peakmode") {
          PARAM_NPA[x0009, 2] <- parallelizationMode
        } else {
          FSA_message("ERROR!!! Problem with NPA0009!")
          checkpoint_parameter <- FALSE
        }
      }
    }
    ##
    RTtolerance <- as.numeric(PARAM_NPA[which(PARAM_NPA[, 1] == 'NPA0010'), 2])
    if (length(RTtolerance) == 0) {
      FSA_message("ERROR!!! Problem with NPA0010! This parameter should be a positive number!")
      checkpoint_parameter <- FALSE
    } else {
      if (RTtolerance <= 0) {
        FSA_message("ERROR!!! Problem with NPA0010! This parameter should be a positive number!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    minSNRbaseline <- as.numeric(PARAM_NPA[which(PARAM_NPA[, 1] == 'NPA0011'), 2])
    if (length(minSNRbaseline) == 0) {
      FSA_message("ERROR!!! Problem with NPA0011! This parameter should be a positive number!")
      checkpoint_parameter <- FALSE
    } else {
      if (minSNRbaseline <= 0) {
        FSA_message("ERROR!!! Problem with NPA0011! This parameter should be a positive number!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    peakHeightThreshold <- as.numeric(PARAM_NPA[which(PARAM_NPA[, 1] == 'NPA0012'), 2])
    if (length(peakHeightThreshold) == 0) {
      FSA_message("ERROR!!! Problem with NPA0012! This parameter should be a positive number!")
      checkpoint_parameter <- FALSE
    } else {
      if (peakHeightThreshold < 0) {
        FSA_message("ERROR!!! Problem with NPA0012! This parameter should be a positive number!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    smoothingWindow <- as.numeric(PARAM_NPA[which(PARAM_NPA[, 1] == 'NPA0013'), 2])
    if (length(smoothingWindow) == 0) {
      FSA_message("ERROR!!! Problem with NPA0013! This parameter should be a positive integer!")
      checkpoint_parameter <- FALSE
    } else {
      if (smoothingWindow <= 0) {
        FSA_message("ERROR!!! Problem with NPA0013! This parameter should be a positive integer!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    nSpline <- as.numeric(PARAM_NPA[which(PARAM_NPA[, 1] == 'NPA0014'), 2])
    if (length(nSpline) == 0) {
      FSA_message("ERROR!!! Problem with NPA0014! This parameter should be a positive integer!")
      checkpoint_parameter <- FALSE
    } else {
      if (nSpline < 10) {
        FSA_message("ERROR!!! Problem with NPA0014! This parameter should be a positive integer!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    topPercentagePeakHeight <- as.numeric(PARAM_NPA[which(PARAM_NPA[, 1] == 'NPA0015'), 2])
    if (length(topPercentagePeakHeight) == 0) {
      FSA_message("ERROR!!! Problem with NPA0015! This parameter should be a positive numberbetween 50-100!")
      checkpoint_parameter <- FALSE
    } else {
      if (topPercentagePeakHeight < 50 | topPercentagePeakHeight > 100) {
        FSA_message("ERROR!!! Problem with NPA0015! This parameter should be a positive number between 50-100!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    minIonRangeDifference <- as.numeric(PARAM_NPA[which(PARAM_NPA[, 1] == 'NPA0016'), 2])
    if (length(minIonRangeDifference) == 0) {
      FSA_message("ERROR!!! Problem with NPA0016! This parameter should be a positive number greater than 0!")
      checkpoint_parameter <- FALSE
    } else {
      if (minIonRangeDifference < 0) {
        FSA_message("ERROR!!! Problem with NPA0016! This parameter should be a positive number greater than 0!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    minNumNPApeaks <- as.numeric(PARAM_NPA[which(PARAM_NPA[, 1] == 'NPA0017'), 2])
    if (length(minNumNPApeaks) == 0) {
      FSA_message("ERROR!!! Problem with NPA0017! This parameter should be a positive number >= 2!")
      checkpoint_parameter <- FALSE
    } else {
      minNumNPApeaks <- ceiling(minNumNPApeaks)
      if (minNumNPApeaks < 2) {
        FSA_message("ERROR!!! Problem with NPA0017! This parameter should be a positive number >= 2!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    pearsonRHOthreshold <- as.numeric(PARAM_NPA[which(PARAM_NPA[, 1] == 'NPA0018'), 2])
    if (length(pearsonRHOthreshold) == 0) {
      FSA_message("ERROR!!! Problem with NPA0018! This parameter should be a positive integer between 0.5-1!")
      checkpoint_parameter <- FALSE
    } else {
      if (pearsonRHOthreshold < 0.5 | pearsonRHOthreshold > 1) {
        FSA_message("ERROR!!! Problem with NPA0018! This parameter should be a positive integer between 0.5-1!")
        checkpoint_parameter <- FALSE
      }
    }
    ############################################################################
    #### Unique tag aggregation by spectra similarity across entire samples ####
    ############################################################################
    x0020 <- which(PARAM_NPA[, 1] == 'NPA0020')
    if (length(x0020) == 0) {
      FSA_message("ERROR!!! Problem with NPA0020!")
      checkpoint_parameter <- FALSE
    } else {
      NPA0020 <- tolower(PARAM_NPA[x0020, 2])
      if (NPA0020 == "yes" | NPA0020 == "no") {
        PARAM_NPA[x0020, 2] <- NPA0020
      } else {
        FSA_message("ERROR!!! Problem with NPA0020!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    if (refMSPcreationCheck) {
      RTtoleranceRef <- as.numeric(PARAM_NPA[which(PARAM_NPA[, 1] == 'NPA0021'), 2])
      if (length(RTtoleranceRef) == 0) {
        FSA_message("ERROR!!! Problem with NPA0021! This parameter should be a positive number!")
        checkpoint_parameter <- FALSE
      } else {
        if (RTtoleranceRef < 0) {
          FSA_message("ERROR!!! Problem with NPA0021! This parameter should be a positive number!")
          checkpoint_parameter <- FALSE
        }
      }
    } else {
      minNPAdetectionFrequency <- as.numeric(PARAM_NPA[which(PARAM_NPA[, 1] == 'NPA0022'), 2])
      if (length(minNPAdetectionFrequency) == 0) {
        FSA_message("ERROR!!! Problem with NPA0022! This parameter should be a positive number between 0 - 100!")
        checkpoint_parameter <- FALSE
      } else {
        if (!((minNPAdetectionFrequency > 0) & (minNPAdetectionFrequency < 100))) {
          FSA_message("ERROR!!! Problem with NPA0022! This parameter should be a positive number between 0 - 100!")
          checkpoint_parameter <- FALSE
        }
      }
    }
    ##
    x0023 <- which(PARAM_NPA[, 1] == 'NPA0023')
    allowedWeightedSpectralEntropy <- tolower(gsub(" ", "", PARAM_NPA[x0023, 2]))
    if (allowedWeightedSpectralEntropy == "1" | allowedWeightedSpectralEntropy == "t" | allowedWeightedSpectralEntropy == "true") {
      allowedWeightedSpectralEntropy <- TRUE
    } else {
      allowedWeightedSpectralEntropy <- FALSE
    }
    PARAM_NPA[x0023, 2] <- allowedWeightedSpectralEntropy
    ##
    minEntropySimilarity <- as.numeric(PARAM_NPA[which(PARAM_NPA[, 1] == 'NPA0024'), 2])
    if (length(minEntropySimilarity) == 0) {
      FSA_message("ERROR!!! Problem with NPA0024! This parameter should be a positive number between 0 - 1!")
      checkpoint_parameter <- FALSE
    } else {
      if (!((minEntropySimilarity >= 0) & (minEntropySimilarity <= 1))) {
        FSA_message("ERROR!!! Problem with NPA0024! This parameter should be a positive number between 0 - 1!")
        checkpoint_parameter <- FALSE
      }
    }
  }
  ##############################################################################
  if (!checkpoint_parameter) {
    PARAM_NPA <- NULL
  }
  return(PARAM_NPA)
}
