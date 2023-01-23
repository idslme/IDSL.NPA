IDSL.NPA_xlsxAnalyzer <- function(spreadsheet) {
  FSA_message("Initiated testing the IDSL.NPA workflow spreadsheet consistency!", failedMessage = FALSE)
  ##
  checkpoint_parameter <- FALSE
  ##
  if (length(spreadsheet) >= 4) {
    if (typeof(spreadsheet) == "list") {
      PARAM_Start <- cbind(spreadsheet[, 2], spreadsheet[, 4])
      checkpoint_parameter <- TRUE
    } else {
      FSA_message("The IDSL.NPA workflow spreadsheet tab was not produced properly!")
      checkpoint_parameter <- FALSE
    }
  } else if (length(spreadsheet) == 1) {
    if (typeof(spreadsheet) == "character") {
      if (file.exists(spreadsheet)) {
        PARAM_Start <- readxl::read_xlsx(spreadsheet, sheet = "Start")
        PARAM_Start <- cbind(PARAM_Start[, 2], PARAM_Start[, 4])
        checkpoint_parameter <- TRUE
      } else {
        FSA_message("The IDSL.NPA workflow spreadsheet not found! It should be an Excel file with .xlsx extention!")
      }
    } else {
      FSA_message("The IDSL.NPA workflow spreadsheet was not produced properly!")
    }
  } else {
    FSA_message("The IDSL.NPA workflow spreadsheet was not produced properly!")
  }
  ##
  if (checkpoint_parameter) {
    ############################################################################
    x0001 <- which(PARAM_Start[, 1] == 'PARAM0001')
    if (length(x0001) == 0) {
      FSA_message("ERROR!!! Problem with PARAM0001!")
      checkpoint_parameter <- FALSE
    } else {
      PARAM0001 <- tolower(PARAM_Start[x0001, 2])
      if (PARAM0001 == "yes" | PARAM0001 == "no") {
        PARAM_Start[x0001, 2] <- PARAM0001
      } else {
        FSA_message("ERROR!!! Problem with PARAM0001!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    x0002 <- which(PARAM_Start[, 1] == 'PARAM0002')
    if (length(x0002) == 0) {
      FSA_message("ERROR!!! Problem with PARAM0002!")
      checkpoint_parameter <- FALSE
    } else {
      PARAM0002 <- tolower(PARAM_Start[x0002, 2])
      if (PARAM0002 == "yes" | PARAM0002 == "no") {
        PARAM_Start[x0002, 2] <- PARAM0002
      } else {
        FSA_message("ERROR!!! Problem with PARAM0002!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    x0003 <- which(PARAM_Start[, 1] == 'PARAM0003')
    if (length(x0003) == 0) {
      FSA_message("ERROR!!! Problem with PARAM0003!")
      checkpoint_parameter <- FALSE
    } else {
      PARAM0003 <- tolower(PARAM_Start[x0003, 2])
      if (PARAM0003 == "yes" | PARAM0003 == "no") {
        PARAM_Start[x0003, 2] <- PARAM0003
      } else {
        FSA_message("ERROR!!! Problem with PARAM0003!")
        checkpoint_parameter <- FALSE
      }
    }
    ############################################################################
    if (!checkpoint_parameter) {
      FSA_message("ERROR!!! Problem with the `Start` spreadsheet tab!")
    }
  }
  ##############################################################################
  ##############################################################################
  ##############################################################################
  if (checkpoint_parameter) {
    ##
    if (PARAM0001 == "yes") {
      FSA_message("Initiated testing the `NPA` spreadsheet tab consistency!", failedMessage = FALSE)
      ##
      PARAM_NPA <- NPA_xlsxAnalyzer(spreadsheet)
      if (is.null(PARAM_NPA)) {
        FSA_message("ERROR!!! Problem with the `NPA` spreadsheet tab!")
        checkpoint_parameter <- FALSE
      }
    } else {
      PARAM_NPA <- NULL
    }
    ##
    if (PARAM0002 == "yes") {
      FSA_message("Initiated testing the `FSDB` spreadsheet tab consistency!", failedMessage = FALSE)
      ##
      PARAM_FSdb <- readxl::read_xlsx(spreadsheet, sheet = "FSDB")
      PARAM_FSdb <- cbind(PARAM_FSdb[, 2], PARAM_FSdb[, 4])
      PARAM_FSdb <- rbind(PARAM_FSdb, c('FSdb0006', TRUE))
      PARAM_FSdb <- FSA_FSdb_xlsxAnalyzer(PARAM_FSdb)
      if (is.null(PARAM_FSdb)) {
        FSA_message("ERROR!!! Problem with the `FSDB` spreadsheet tab!")
        checkpoint_parameter <- FALSE
      }
    } else {
      PARAM_FSdb <- NULL
    }
    ##
    if (PARAM0003 == "yes") {
      FSA_message("Initiated testing the `SpectraSimilarity` spreadsheet tab consistency!", failedMessage = FALSE)
      ##
      PARAM_SPEC <- readxl::read_xlsx(spreadsheet, sheet = "SpectraSimilarity")
      PARAM_SPEC <- cbind(PARAM_SPEC[, 2], PARAM_SPEC[, 4])
      PARAM_SPEC <- rbind(PARAM_SPEC, c('SPEC0003', TRUE))
      PARAM_SPEC <- FSA_SpectraSimilarity_xlsxAnalyzer(PARAM_SPEC)
      if (is.null(PARAM_SPEC)) {
        FSA_message("ERROR!!! Problem with the `SpectraSimilarity` spreadsheet tab!")
        checkpoint_parameter <- FALSE
      }
    } else {
      PARAM_SPEC <- NULL
    }
    #################### MS/MS library import and export #######################
    x0004 <- which(PARAM_Start[, 1] == 'PARAM0004')
    ##
    if (length(x0004) == 0) {
      FSA_message("ERROR!!! Problem with PARAM0004!")
      checkpoint_parameter <- FALSE
    } else {
      ##
      if ((PARAM0002 == "no") & (PARAM0003 == "yes")) {
        ##
        FSdb_file <- PARAM_Start[x0004, 2]
        FSdb_file <- gsub("\\", "/", FSdb_file, fixed = TRUE)
        PARAM_Start[x0004, 2] <- FSdb_file
        if (!file.exists(FSdb_file)) {
          FSA_message("ERROR!!! Problem with PARAM0004! Please ensure the full path is provided for the FSDB in .Rdata format OR select 'YES' for PARAM0002!")
          checkpoint_parameter <- FALSE
        }
      }
    }
    ############################################################################
    ####################### Sample import and export ###########################
    if (PARAM0003 == "yes") {
      ##
      x0005 <- which(PARAM_Start[, 1] == 'PARAM0005')
      x0006 <- which(PARAM_Start[, 1] == 'PARAM0006')
      ##
      if (PARAM0001 == "no") {
        if (length(x0005) == 0) {
          FSA_message("ERROR!!! Problem with PARAM0005!")
          checkpoint_parameter <- FALSE
        } else {
          address_input_msp <- PARAM_Start[x0005, 2]
          address_input_msp <- gsub("\\", "/", address_input_msp, fixed = TRUE)
          PARAM_Start[x0005, 2] <- address_input_msp
          if (!dir.exists(address_input_msp)) {
            FSA_message("ERROR!!! Problem with PARAM0005! Please make sure the full path is provided!")
            checkpoint_parameter <- FALSE
          } else {
            file_name_sample_msp <- dir(path = address_input_msp)
            file_name_sample_msp <- file_name_sample_msp[grepl(".msp$", file_name_sample_msp, ignore.case = TRUE)]
            if (length(file_name_sample_msp) == 0) {
              FSA_message("ERROR!!! Problem with PARAM0005! No .msp file was detected in the designated folder!")
              checkpoint_parameter <- FALSE
            }
          }
        }
        ########################################################################
        if (length(x0006) == 0) {
          FSA_message("ERROR!!! Problem with PARAM0006!")
          checkpoint_parameter <- FALSE
        } else {
          output_sample <- PARAM_Start[x0006, 2]
          output_sample <- gsub("\\", "/", output_sample, fixed = TRUE)
          PARAM_Start[x0006, 2] <- output_sample
          if (!dir.exists(output_sample)) {
            tryCatch(dir.create(output_sample, recursive = TRUE), warning = function(w){warning("Problem with PARAM0006! R cannot create the folder!")})
            if (!dir.exists(output_sample)) {
              checkpoint_parameter <- FALSE
            }
          }
        }
        ##
      } else if (PARAM0001 == "yes") {
        if (!is.null(PARAM_NPA)) {
          ##
          refMSPcreationCheck <- file.exists(as.character(PARAM_NPA[which(PARAM_NPA[, 1] == 'NPA0005'), 2]))
          if (refMSPcreationCheck) {
            annex_MSP <- "/CSA_REF_MSP"
          } else {
            annex_MSP <- "/CSA_MSP"
          }
          ##
          xNPA0008 <- which(PARAM_NPA[, 1] == 'NPA0008')
          address_input_msp <- PARAM_NPA[xNPA0008, 2]
          address_input_msp <- gsub("\\", "/", address_input_msp, fixed = TRUE)
          address_input_msp <- paste0(address_input_msp, annex_MSP)
          PARAM_Start[x0005, 2] <- address_input_msp
          ##
          output_path_sample <- PARAM_NPA[xNPA0008, 2]
          output_path_sample <- gsub("\\", "/", output_path_sample, fixed = TRUE)
          PARAM_Start[x0006, 2] <- output_path_sample
        }
      }
    }
    ############################################################################
    if ((PARAM0002 == "yes") & !is.null(PARAM_SPEC) & (PARAM0003 == "yes") & !is.null(PARAM_FSdb)) {
      ##
      FSdb0007 <- as.numeric(PARAM_FSdb[which(PARAM_FSdb[, 1] == 'FSdb0007'), 2])
      SPEC0008 <- as.numeric(PARAM_SPEC[which(PARAM_SPEC[, 1] == 'SPEC0008'), 2])
      if (FSdb0007 != SPEC0008) {
        FSA_message("Inconsistency between `FSdb0007` and `SPEC0008` parameters in the `FSDB` and `SpectraSimilarity` tabs!")
        checkpoint_parameter <- FALSE
      }
      ##
      FSdb0009 <- eval(parse(text = PARAM_FSdb[which(PARAM_FSdb[, 1] == 'FSdb0009'), 2]))
      SPEC0013 <- eval(parse(text = PARAM_SPEC[which(PARAM_SPEC[, 1] == 'SPEC0013'), 2]))
      if (FSdb0009 != SPEC0013) {
        FSA_message("Inconsistency between `FSdb0009` and `SPEC0013` parameters in the `FSDB` and `SpectraSimilarity` tabs!")
        checkpoint_parameter <- FALSE
      }
    }
  }
  ##############################################################################
  if (checkpoint_parameter) {
    ##
    PARAM_total <- list(PARAM_Start, PARAM_NPA, PARAM_FSdb, PARAM_SPEC)
    ##
    if (PARAM0001 == "yes") {
      FSA_message("The `NPA` spreadsheet tab is consistent with the IDSL.NPA workflow!", failedMessage = FALSE)
    }
    ##
    if (PARAM0002 == "yes") {
      FSA_message("The `FSDB` spreadsheet tab is consistent with the IDSL.NPA workflow!", failedMessage = FALSE)
    }
    ##
    if (PARAM0003 == "yes") {
      ##
      massErrorPrecursor <- tryCatch(as.numeric(PARAM_SPEC[which(PARAM_SPEC[, 1] == 'SPEC0005'), 2]), warning = function(w) {NA})
      if (!is.na(massErrorPrecursor)) {
        FSA_message("NOTICE: Precursor m/z match (SPEC0005) was selected for spectra annotations! Empty annotation tables will be generated when .msp files do not contain 'PrecursorMZ' information!", failedMessage = FALSE)
      } else {
        FSA_message("NOTICE: Precursor m/z match (SPEC0005) was not selected for spectra annotations!", failedMessage = FALSE)
      }
      ##
      RTtolerance <- PARAM_SPEC[which(PARAM_SPEC[, 1] == 'SPEC0006'), 2]
      if (!is.na(RTtolerance)) {
        FSA_message(paste0("NOTICE: Retention time tolerance = '" , RTtolerance, " min' (SPEC0006) was selected for spectra annotations! Empty annotation tables will be generated when .msp files do not contain 'Retention Time' information!"), failedMessage = FALSE)
      } else {
        FSA_message("NOTICE: Retention time tolerance (SPEC0006) was not selected and retention time values will not be used for spectra annotations!", failedMessage = FALSE)
      }
      ##
      FSA_message("The `SpectraSimilarity` spreadsheet tab is consistent with the IDSL.NPA workflow!", failedMessage = FALSE)
    }
    ##
    FSA_message("The spreadsheet is consistent with the IDSL.NPA workflow!", failedMessage = FALSE)
  } else {
    ##
    FSA_message("The spreadsheet is not consistent with the IDSL.NPA workflow!")
    PARAM_total <- vector(mode = "list", 4)
  }
  ##
  names(PARAM_total) <- c("PARAM_Start", "PARAM_NPA", "PARAM_FSdb", "PARAM_SPEC")
  ##############################################################################
  return(PARAM_total)
}
