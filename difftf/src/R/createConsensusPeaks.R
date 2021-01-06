start.time  <-  Sys.time()



#########################
# LIBRARY AND FUNCTIONS #
#########################

library("checkmate")
assertClass(snakemake, "Snakemake")
assertDirectoryExists(snakemake@config$par_general$dir_scripts)
source(paste0(snakemake@config$par_general$dir_scripts, "/functions.R"))

initFunctionsScript(packagesReq = NULL, minRVersion = "3.1.0", warningsLevel = 1, disableScientificNotation = TRUE)
checkAndLoadPackages(c("tidyverse", "futile.logger", "DiffBind", "checkmate", "stats"), verbose = FALSE)

########################################################################
# SAVE SNAKEMAKE S4 OBJECT THAT IS PASSED ALONG FOR DEBUGGING PURPOSES #
########################################################################

# Use the following line to load the Snakemake object to manually rerun this script (e.g., for debugging purposes)
# Replace {outputFolder} correspondingly.
# snakemake = readRDS("{outputFolder}/LOGS_AND_BENCHMARKS/createConsensusPeaks.R.rds")
createDebugFile(snakemake)

###################
#### PARAMETERS ###
###################

par.l = list()

par.l$verbose = TRUE
par.l$log_minlevel = "INFO"

#####################
# VERIFY PARAMETERS #
#####################


assertClass(snakemake, "Snakemake")

## INPUT ##
assertList(snakemake@input, min.len = 1)
assertSubset(names(snakemake@input), c("", "peaks", "checkFlag", "sampleFile"))

par.l$file_input_peaks = snakemake@input$peaks


## OUTPUT ##
assertList(snakemake@output, min.len = 1)
assertSubset(names(snakemake@output), c("", "consensusPeaks_bed", "summaryPlot"))

par.l$file_output_consensusPeaks = snakemake@output$consensusPeaks_bed
par.l$file_output_plot           = snakemake@output$summaryPlot


## CONFIG ##
assertList(snakemake@config, min.len = 1)

file_sampleData = snakemake@config$samples$summaryFile
assertFileExists(file_sampleData)

peakType = snakemake@config$peaks$peakType
assertSubset(peakType, c("raw", "bed", "narrow", "macs", "swembl", "bayes", "peakset", "fp4"), empty.ok = FALSE)

minOverlap = snakemake@config$peaks$minOverlap
assertIntegerish(minOverlap, lower = 0)

par.l$debugMode = setDebugMode(snakemake@config$par_general$debugMode)


## LOG ##
assertList(snakemake@log, min.len = 1)
par.l$file_log = snakemake@log[[1]]


allDirs = c(dirname(par.l$file_output_consensusPeaks), 
            dirname(par.l$file_output_plot),
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
    sampleMetaData.df = read_tidyverse_wrapper(file_sampleData, type = "tsv", col_types = cols())
  
    save(list = ls(), file = snakemake@params$debugFile)
    flog.info(paste0("File ", snakemake@params$debugFile, " has been saved. You may use it for trouble-shooting and debugging, see the Documentation for more details."))
    
}


##########################
# Create consensus peaks #
##########################

# Provide the metadata file and parse the CSV here
sampleMetaData.df = read_tidyverse_wrapper(file_sampleData, type = "tsv", col_types = cols())
assertSubset(c("SampleID", "bamReads", "conditionSummary", "Peaks"), colnames(sampleMetaData.df))
assertIntegerish(minOverlap, upper = nrow(sampleMetaData.df))

sampleMetaData.df$Condition  = sampleMetaData.df$conditionSummary
sampleMetaData.df$PeakCaller = peakType
sampleMetaData.df$conditionSummary = sampleMetaData.df$Condition

# sampleMetaData.df = sampleMetaData.df[complete.cases(sampleMetaData.df),]

if (nrow(sampleMetaData.df) == 0) {
  message = paste0("No rows left to analyze, the data frame for DiffBind contains no rows after filtering for complete cases")
  checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
}
## DiffBind assigning. Convert to regular data frame first so DiffBind does not complain
dba <- dba(sampleSheet = as.data.frame(sampleMetaData.df))
consensus = dba.peakset(dba, minOverlap = minOverlap)

# Get the consensus peaks (last row)
consensusPeaks.df = as.data.frame(consensus$peaks[[nrow(sampleMetaData.df) + 1]])

colnames(consensusPeaks.df) = c("chr", "start", "end", "pvalue")
consensusPeaks.df$length = consensusPeaks.df$end - consensusPeaks.df$start

# Graph
# TODO: Make a bit nicer
g = ggplot(consensusPeaks.df) + geom_density(aes(x = length)) + scale_x_continuous(labels = seq(0,3000,250), breaks = seq(0,3000,250))
ggsave(plot = g, filename = par.l$file_output_plot)


consensusPeaks.df$annotation = paste0(consensusPeaks.df$chr, ":", consensusPeaks.df$start, "-", consensusPeaks.df$end)

# Deactivated now
# consensusPeaks.df$strand = "+"
# Reorganize columns
consensusPeaks.df = consensusPeaks.df[,c("chr", "start", "end", "annotation", "pvalue")]

consensusPeaks.df = mutate_if(consensusPeaks.df, is.numeric, as.character)
write_tsv(consensusPeaks.df, path = par.l$file_output_consensusPeaks, col_names = FALSE)

.printExecutionTime(start.time)

flog.info("Session info: ", sessionInfo(), capture = TRUE)
