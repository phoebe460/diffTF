start.time  <-  Sys.time()


#########################
# LIBRARY AND FUNCTIONS #
#########################
library("checkmate")
assertClass(snakemake, "Snakemake")
assertDirectoryExists(snakemake@config$par_general$dir_scripts)
source(paste0(snakemake@config$par_general$dir_scripts, "/functions.R"))

initFunctionsScript(packagesReq = NULL, minRVersion = "3.1.0", warningsLevel = 1, disableScientificNotation = TRUE)
checkAndLoadPackages(c("tidyverse", "futile.logger", "modeest", "checkmate", "ggrepel"), verbose = FALSE)

########################################################################
# SAVE SNAKEMAKE S4 OBJECT THAT IS PASSED ALONG FOR DEBUGGING PURPOSES #
########################################################################

# Use the following line to load the Snakemake object to manually rerun this script (e.g., for debugging purposes)
# Replace {outputFolder} correspondingly.
# snakemake = readRDS("{outputFolder}/LOGS_AND_BENCHMARKS/summary1.R.rds")
createDebugFile(snakemake)

###################
#### PARAMETERS ###
###################
par.l = list()

# Hard-coded parameters
par.l$verbose = TRUE
par.l$min_pValue   = .Machine$double.xmin
par.l$FDR_threshold = 0.05
par.l$plot_min_diffMean = 0.02
par.l$log_minlevel = "INFO"


#####################
# VERIFY PARAMETERS #
#####################

assertClass(snakemake, "Snakemake")

## INPUT ##
assertList(snakemake@input, min.len = 1)
assertSubset(names(snakemake@input), c("", "peaks", "TF", "analysesTFOutput"))

par.l$file_input_peaks = snakemake@input$peaks
assertFileExists(par.l$file_input_peaks, access = "r")

par.l$files_input_TF_summary = snakemake@input$TF
for (fileCur in par.l$files_input_TF_summary) {
  assertFileExists(fileCur, access = "r")
}

## OUTPUT ##
assertList(snakemake@output, min.len = 1)
assertSubset(names(snakemake@output), c("", "outputTable"))

#par.l$file_output_volcanoPlot = snakemake@output$volcanoPlot
par.l$file_output_table       = snakemake@output$outputTable

par.l$debugMode = setDebugMode(snakemake@config$par_general$debugMode)

## LOG ##
assertList(snakemake@log, min.len = 1)
par.l$file_log = snakemake@log[[1]]



allDirs = c(dirname(par.l$file_output_table),
            dirname(par.l$file_log)
          )

testExistanceAndCreateDirectoriesRecursively(allDirs)


######################
# FINAL PREPARATIONS #
######################

startLogger(par.l$file_log, par.l$log_minlevel, removeOldLog = TRUE)
printParametersLog(par.l)

if (par.l$debugMode) {
    
    flog.info(paste0("Debug mode active. Reading it all files that this step requires and save it to ", snakemake@params$debugFile))
    
    
    # Read all files already here and then save the session so as much as possible from the script can be executed without file dependencies
    peaks.df = read_tidyverse_wrapper(par.l$file_input_peaks, type = "tsv", col_types = cols())

    par.l$filesTransl.l = list()
    nTF = length(par.l$files_input_TF_summary)
    
    # Create translation table
    for (fileCur in par.l$files_input_TF_summary) {
        elements = strsplit(fileCur, split = "/",  fixed = TRUE)[[1]]
        hit =  which(grepl(pattern = "^extension", elements))
        stopifnot(length(hit) == 1 & hit != 1)
        par.l$filesTransl.l[[fileCur]] = elements[hit - 1]  
    }
    
    stats.l = list()
    for (fileCur in par.l$files_input_TF_summary) {
        
        TF = par.l$filesTransl.l[[fileCur]]
        stopifnot(!is.null(TF))
        if (file.exists(fileCur)) {
            stats.l[[fileCur]] = readRDS(fileCur)
        } else {
            message = paste0("File missing: ", fileCur)
            checkAndLogWarningsAndErrors(NULL,  message, isWarning = TRUE)
        }
    }
    
    save(list = ls(), file = snakemake@params$debugFile)
    flog.info(paste0("File ", snakemake@params$debugFile, " has been saved. You may use it for trouble-shooting and debugging, see the Documentation for more details."))
    
}



################
# Collect data #
################

par.l$filesTransl.l = list()
nTF = length(par.l$files_input_TF_summary)

# Create translation table
for (fileCur in par.l$files_input_TF_summary) {
  elements = strsplit(fileCur, split = "/",  fixed = TRUE)[[1]]
  hit =  which(grepl(pattern = "^extension", elements))
  stopifnot(length(hit) == 1 & hit != 1)
  par.l$filesTransl.l[[fileCur]] = elements[hit - 1]  
}


peaks.df = read_tidyverse_wrapper(par.l$file_input_peaks, type = "tsv", col_types = cols())

summary.df = NULL

nTFs = length(par.l$filesTransl.l)

flog.info(paste0("Using ", nTFs, " TFs"))

for (fileCur in par.l$files_input_TF_summary) {
  
  TF = par.l$filesTransl.l[[fileCur]]
  stopifnot(!is.null(TF))
  if (file.exists(fileCur)) {
    
    stats.df = readRDS(fileCur)
    
    if (nrow(stats.df) == 0) {
        message = paste0("File ", fileCur, " contains no rows, this TF will be skipped thereafter")
        checkAndLogWarningsAndErrors(NULL,  message, isWarning = TRUE)
    }
    
    if (is.null(summary.df)) {
      summary.df = stats.df
    } else {
      summary.df = rbind(summary.df, stats.df)
    }
    

  } else {
    message = paste0("File missing: ", fileCur)
    checkAndLogWarningsAndErrors(NULL,  message, isWarning = TRUE)
  }
  
}

nTF = length(unique(summary.df$TF))

flog.info(paste0(" Imported ", nTF, " TFs out of a list of ", nTFs, ". Missing: ", nTFs - nTF))

nTFMissing = length(which(is.na(summary.df$Pos_l2FC)))
if (nTFMissing == nrow(summary.df)) {
  error = "All TF have missing data. Cannot continue. Add more samples or change the peaks."
  checkAndLogWarningsAndErrors(NULL,  error, isWarning = FALSE)
}
  

# Replace p-values of 0 with the smallest p-value on the system
summary.df$pvalue_raw[summary.df$pvalue_raw == 0] = .Machine$double.xmin

# Check the version of modeest, because version 2.3.2 introduced an implementation change that breaks things

if (packageVersion("modeest") < "2.3.2") {
    
    mode_peaks = mlv(peaks.df$l2FC, method = "mfv", na.rm = TRUE)
    
    stopifnot(is.list(mode_peaks))
    l2fc_mode = ifelse(is.null(mode_peaks$M), NA, mode_peaks$M)
    l2fc_skewness = ifelse(is.null(mode_peaks$skewness), NA, mode_peaks$skewness)
} else {
    
    l2fc_mode = mlv(peaks.df$l2FC, method = "mfv", na.rm = TRUE)[1]
    l2fc_skewness = skewness(peaks.df$l2FC, na.rm = TRUE)[1]
}



summary.df = summary.df %>%
              dplyr::mutate(
                  adj_pvalue  = p.adjust(pvalue_raw, method = "fdr"),
                  Diff_mean   = Mean_l2FC    -   mean(peaks.df$l2FC, na.rm = TRUE), 
                  Diff_median = Median_l2FC  - median(peaks.df$l2FC, na.rm = TRUE),
                  Diff_mode   = Mode_l2FC    -  l2fc_mode,    
                  Diff_skew   = skewness_l2FC -  l2fc_skewness)  %>%
              na.omit(summary.df) %>%
             mutate_at(c("Diff_mean", "Diff_median", "Diff_mode", "Diff_skew", "Mean_l2FC", 
                        "Median_l2FC", "Mode_l2FC", "skewness_l2FC"), signif, 3) %>%
             mutate_at(c("pvalue_raw", "adj_pvalue"), formatC, format = "g", digits = 3)


# summary.df = mutate_if(summary.df, is.numeric, as.character)
write_tsv(summary.df, par.l$file_output_table) # TODO: check the dec = "." parameter

.printExecutionTime(start.time)

flog.info("Session info: ", sessionInfo(), capture = TRUE)