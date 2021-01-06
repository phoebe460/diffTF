start.time  <-  Sys.time()

#########################
# LIBRARY AND FUNCTIONS #
#########################

library("checkmate")
assertClass(snakemake, "Snakemake")
assertDirectoryExists(snakemake@config$par_general$dir_scripts)
source(paste0(snakemake@config$par_general$dir_scripts, "/functions.R"))

initFunctionsScript(packagesReq = NULL, minRVersion = "3.1.0", warningsLevel = 1, disableScientificNotation = TRUE)
checkAndLoadPackages(c("tidyverse", "futile.logger", "checkmate", "Rsamtools"), verbose = FALSE)

########################################################################
# SAVE SNAKEMAKE S4 OBJECT THAT IS PASSED ALONG FOR DEBUGGING PURPOSES #
########################################################################

# Use the following line to load the Snakemake object to manually rerun this script (e.g., for debugging purposes)
# Replace {outputFolder} correspondingly.
# snakemake = readRDS("{outputFolder}/LOGS_AND_BENCHMARKS/checkParameters.R.rds")
createDebugFile(snakemake)

par.l = list()

par.l$verbose = TRUE
par.l$log_minlevel = "INFO"


#####################
# VERIFY PARAMETERS #
#####################

assertClass(snakemake, "Snakemake")

## OUTPUT ##
assertList(snakemake@output, min.len = 1)
assertSubset(names(snakemake@output), c("", "consPeaks", "flag"))

par.l$output_peaksClean = snakemake@output$consPeaks
par.l$output_flag = snakemake@output$flag

## CONFIG ##
assertList(snakemake@config, min.len = 1)

file_peaks = snakemake@config$peaks$consensusPeaks

par.l$designFormula = snakemake@config$par_general$designContrast
checkAndLogWarningsAndErrors(par.l$designFormula, checkCharacter(par.l$designFormula, len = 1, min.chars = 3))


TFBS_dir = snakemake@config$additionalInputFiles$dir_TFBS
assertDirectoryExists(dirname(TFBS_dir), access = "r")

fastaFile = snakemake@config$additionalInputFiles$refGenome_fasta
assertFileExists(fastaFile)
assertDirectoryExists(dirname(fastaFile), access = "r")

allTFs = strsplit(snakemake@config$par_general$TFs, ",")[[1]]


par.l$nPermutations = snakemake@config$par_general$nPermutations
assertIntegerish(par.l$nPermutations, lower = 0, len = 1)

par.l$nBootstraps = as.integer(snakemake@config$par_general$nBootstraps)
assertIntegerish(par.l$nBootstraps, len = 1)

par.l$file_input_sampleData = snakemake@config$samples$summaryFile
checkAndLogWarningsAndErrors(par.l$file_input_sampleData, checkFileExists(par.l$file_input_sampleData, access = "r"))

par.l$conditionComparison  = snakemake@config$par_general$conditionComparison
checkAndLogWarningsAndErrors(par.l$conditionComparison, checkCharacter(par.l$conditionComparison, len = 1))

par.l$plotRNASeqClassification = as.logical(snakemake@config$par_general$RNASeqIntegration)
assertFlag(par.l$plotRNASeqClassification)

par.l$debugMode = setDebugMode(snakemake@config$par_general$debugMode)

### PARAMS ##
par.l$TFBSPattern = snakemake@params$suffixTFBS

## LOG ##
assertList(snakemake@log, min.len = 1)
par.l$file_log = snakemake@log[[1]]


allDirs = c(dirname(par.l$file_log))

testExistanceAndCreateDirectoriesRecursively(allDirs)


######################
# FINAL PREPARATIONS #
######################
startLogger(par.l$file_log, par.l$log_minlevel, removeOldLog = TRUE)
printParametersLog(par.l)

if (par.l$debugMode) {
    
    flog.info(paste0("Debug mode active. Reading it all files that this step requires and save it to ", snakemake@params$debugFile))
    
    
    # Read all files already here and then save the session so as much as possible from the script can be executed without file dependencies
    sampleData.df = read_tidyverse_wrapper(par.l$file_input_sampleData, type = "tsv", col_names = TRUE, col_types = cols())
    
    if (file_peaks != "") {
        peaks.df = read_tidyverse_wrapper(snakemake@config$peaks$consensusPeaks, type = "tsv", col_names = FALSE)
    }
   
    
    
    useAllTFs = FALSE 
    if (length(allTFs) == 1) 
        if (allTFs == "all") {
            useAllTFs = TRUE
        }
    
    if (useAllTFs) {
        
        TFs = createFileList(TFBS_dir, par.l$TFBSPattern, verbose = FALSE)
        
        if (length(TFs) == 0) {
            message = paste0("No files with pattern ", par.l$TFBSPattern, " found in directory ", TFBS_dir, " as specified by the parameter \"dir_TFBS\"")
            checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
        }
        
        allTFs = gsub(par.l$TFBSPattern, "", basename(TFs))
        stopifnot(length(allTFs) > 0)
        
    }
    
    TFs.l = list()
    for (TFCur in allTFs) {
        
        TFCur = gsub(pattern = " ",replacement = "",TFCur)
        flog.info(paste0(" Checking TF ", TFCur, "..."))
        fileCur = paste0(TFBS_dir, "/", TFCur, par.l$TFBSPattern)
        if (!file.exists(fileCur)) {
            message = paste0("File ", fileCur, " does not exist even though the TF ", TFCur, " has been specified")
            checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
        }
        
        TFs.l[[TFCur]] = read_tidyverse_wrapper(fileCur, type = "tsv", ncolExpected = 6, col_names = FALSE, col_types = "ciicnc", minRows = 1)
 
    }
    
    save(list = ls(), file = snakemake@params$debugFile)
    flog.info(paste0("File ", snakemake@params$debugFile, " has been saved. You may use it for trouble-shooting and debugging, see the Documentation for more details."))
    
}



######################
# LOAD ALL LIBRARIES #
######################

# Step 1
# TODO: First loading DESeq2 before DiffBind seems to prevent the segfault
checkAndLoadPackages(c("tidyverse", "futile.logger", "DESeq2", "DiffBind", "checkmate", "stats"), verbose = FALSE)

# Step 2: diffPeaks
checkAndLoadPackages(c("tidyverse", "futile.logger", "DESeq2", "vsn", "csaw", "checkmate", "limma", "tools", "geneplotter", "RColorBrewer", "matrixStats"), verbose = FALSE)

# Step 3: analyzeTF
checkAndLoadPackages(c("tidyverse", "futile.logger", "DESeq2", "vsn", "modeest", "checkmate", "limma", "geneplotter", "RColorBrewer", "tools"), verbose = FALSE)

# Step 4: summary1
checkAndLoadPackages(c("tidyverse", "futile.logger", "modeest", "checkmate", "ggrepel"), verbose = FALSE)

# Step 5: binningTF
checkAndLoadPackages(c("tidyverse", "futile.logger", "checkmate", "tools", "methods", "boot", "lsr"), verbose = FALSE)

# Step 6: summaryFinal
checkAndLoadPackages(c("tidyverse", "futile.logger","ggrepel", "checkmate", "tools", "methods", "grDevices", "pheatmap"), verbose = FALSE)

# Only needed if RNA-Seq integration is activated
if (par.l$plotRNASeqClassification) {
    checkAndLoadPackages(c( "lsr", "DESeq2",  "matrixStats",  "pheatmap", "preprocessCore", "apeglm"), verbose = FALSE)
}

# Check the version of readr, at least 1.1.0 is required to properly write gz files

if (packageVersion("readr") < "1.1.0") {
  message = paste0("Version of readr library is too old, at least version 1.1.0 is required). Execute the following in R: install.packages('readr') ") 
  checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
}


if (par.l$nPermutations == 0 && par.l$nBootstraps < 1000) {
  flog.warn(paste0("The value for nBootstraps is < 1000. We strongly recommend using a higher value in order to reliably estimate the statistical variance."))
}

#############################
# CHECK FASTA AND BAM FILES #
#############################

sampleData.df = read_tidyverse_wrapper(par.l$file_input_sampleData, type = "tsv", col_names = TRUE, col_types = cols())

components3types = checkDesignIntegrity(snakemake, par.l, sampleData.df)$types


# Check the sample table. Distinguish between the two modes: Factor and integer for conditionSummary
if (components3types["conditionSummary"] == "logical" | components3types["conditionSummary"] == "factor") {
    
    conditionComparison = strsplit(snakemake@config$par_general$conditionComparison, ",")[[1]]
    assertVector(conditionComparison, len = 2)
    
    nDistValues = length(unique(sampleData.df$conditionSummary))
    if (nDistValues != 2) {
        message = paste0("The column 'conditionSummary' must contain exactly 2 different values, but ", nDistValues, " were found.") 
        checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
    }
    
    if (!testSubset(unique(sampleData.df$conditionSummary), conditionComparison)) {
        message = paste0("The elements specified in 'conditionComparison' in the config file must be a subset of the values in the column 'conditionSummary' in the sample file") 
        checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
    }
} else {
    
    # Test whether at least two distinct values are present
    nDistinct = length(unique(sampleData.df$conditionSummary))
    if (nDistinct < 2) {
        message = paste0("At least 2 distinct values for the column 'conditionSummary' must be present in the sample file in the quantitative mode.") 
        checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
    }
}


# Check the design formula

designFormula = convertToFormula(par.l$designFormula, colnames(sampleData.df))
designMatrix = model.matrix(designFormula, data = sampleData.df)

if (nrow(designMatrix) < nrow(sampleData.df)) {
  missingRows = setdiff(1:nrow(sampleData.df), as.integer(row.names(designMatrix)))
  message = paste0("There is a problem with the specified design formula (parameter designContrast): The corresponding design matrix has fewer rows. This usually means that there are missing values in one of the specified variables. The problem comes from the following lines in the summary file: ", paste0(missingRows, collapse = ","), ".") 
  checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
}

dummyMatrix = matrix(nrow = 100, ncol = nrow(sampleData.df), sample(1:100, 100 * nrow(sampleData.df), replace = TRUE))

TF.DESeq.obj <- DESeqDataSetFromMatrix(countData = dummyMatrix,
                                       colData = sampleData.df,
                                       design = designFormula)

# Build the fasta index. Requires write access to the folder where the fasta is stored (limitation of samtools faidx)

fastaIndex = paste0(fastaFile, ".fai")
if (!file.exists(fastaIndex)) {
  flog.info(paste0("Running samtools faidx for fasta file to generate fasta index"))
  indexFa(fastaFile)
} else {
  flog.info(paste0("Fasta index already found"))
}

indexes.df = as.data.frame(scanFaIndex(fastaFile))
indexes.df$seqnames = as.character(indexes.df$seqnames)

# Check individually for each BAM file
for (bamCur in sampleData.df$bamReads) {
  
  bamHeader = scanBamHeader(bamCur)
  
  chrBAM.df = as.data.frame(bamHeader[[1]]$targets)
  chrBAM.df$seqnames = rownames(chrBAM.df)
  colnames(chrBAM.df) = c("width", "seqnames")
  
  # Check whether BAM file contains chromosomes with the "chr" notation
  if (length(which(grepl("^chr", chrBAM.df$seqnames))) == 0) {
      message = paste0("File ", bamCur, " does not have the correct chromosome names. The \"chr\" prefix is required for proper chromosome names, but they were not found. Check your BAM files and use samtools and sed to add \"chr\" to each chromosome name") 
      checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
  }
  
  merged.df = full_join(chrBAM.df, indexes.df, by = "seqnames")
  
  # Filter chromosomes and retain only those we keep in the pipeline
  discardMatches = grepl("^chrX|^chrY|^chrM|^chrUn|random|hap|_gl|^GL", merged.df$seqnames, perl = TRUE)
  
  if (length(discardMatches) > 0) {
    chrNamesDiscarded = paste0(merged.df$seqnames[discardMatches], collapse = " ", sep = "\n")
    flog.info(paste0("Discard the following chromosomes for compatibility comparisons because they are filtered in later steps anyway:\n ", chrNamesDiscarded))
    merged.df =  merged.df[!discardMatches,]
  }
  
  
  mismatches = merged.df$seqnames[which(merged.df$width.x != merged.df$width.y)]
  if (length(mismatches) > 0) {
    chrStr = paste0(mismatches, collapse = ",")
    message = paste0("Chromosome lengths differ between the fasta file ", fastaFile, " and the BAM file ", bamCur, " for chromosome(s) ", chrStr) 
    checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
  } else {
    flog.info(paste0("Reference fasta file ", fastaFile, " and BAM file ", bamCur, " are fully compatible."))
  }
  
}

flog.info(paste0("Check peak files..."))

###################
# CHECK PEAK FILE #
###################

if (file_peaks != "") {
  
  assertFileExists(snakemake@config$peaks$consensusPeaks)
  peaks.df = read_tidyverse_wrapper(snakemake@config$peaks$consensusPeaks, type = "tsv", col_names = FALSE)

  flog.info(paste0("Peak file contains ", nrow(peaks.df), " peaks."))
  
  if (nrow(peaks.df) > 100000) {
    message = paste0("The number of peaks is high (", nrow(peaks.df), "), subsequent steps may be slow, particularly the binningTF steps. If you encounter problems related to execution time during the analysis, consider decreasing the number of peaks for improved performance.")
    checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
    
    if (snakemake@config$par_general$nPermutations > 5) {
      message = paste0("In addition to the high number of peaks, more than 5 permutations have been selected, which will further increase running times and memory footprint. Consider decreasing the number of peaks for improved performance.")
      checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
    }
  }
  
  if (ncol(peaks.df) < 3) {
    message = paste0("At least 3 columns are required, but only ", ncol(peaks.df), " columns have been found in the file ", snakemake@config$peaks$consensusPeaks)
    checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
  }
  
  assertCharacter(peaks.df$X1)
  assertIntegerish(peaks.df$X2)
  assertIntegerish(peaks.df$X3)
  
  # End coordinates must not be smaller than start coordinates
  assertIntegerish(peaks.df$X3 - peaks.df$X2, lower = 0)
  
  if (ncol(peaks.df) == 3) {
    peaks.df$X4 = "TODO"
    peaks.df$X5 = 0
  } else if (ncol(peaks.df) == 4) {
    peaks.df$X5 = 0
  } else if (ncol(peaks.df) == 5) {
    
  } else {
    message = paste0("More than 5 columns found. Only the first 5 are needed, the remaining ones will be dropped for subsequent analysis...")
    checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
    peaks.df = peaks.df[,1:5]
  }
  
  # Overwrite original annotation, we use our own from here on
  peaks.df$X4 = paste0(peaks.df$X1, ":", peaks.df$X2, "-", peaks.df$X3)
  
  
  colnamesFinal = c("chr", "start", "end", "annotation", "score")
  colnames(peaks.df) = colnamesFinal
  
  rowsBefore = nrow(peaks.df)
  
  peaks.df = dplyr::distinct(peaks.df, annotation, .keep_all = TRUE)
  
  nRowsFiltered = rowsBefore - nrow(peaks.df)
  if (par.l$verbose & nRowsFiltered  > 0) flog.info(paste0("Filtered ", nRowsFiltered, " non-unique positions out of ", rowsBefore, " from peaks."))
  
  peaks.df.transf = dplyr::mutate_if(peaks.df, is.numeric, as.character)
  write_tsv(peaks.df.transf, path = par.l$output_peaksClean, col_names = FALSE)  
  
} else {
  
  cat("DUMMY FILE, DO NOT DELETE", file = par.l$output_peaksClean) 
  
}

flog.info(paste0("Check TF-specific TFBS files..."))

useAllTFs = FALSE 
if (length(allTFs) == 1) 
  if (allTFs == "all") {
    useAllTFs = TRUE
  }

if (useAllTFs) {
  
  TFs = createFileList(TFBS_dir, par.l$TFBSPattern, verbose = FALSE)
  
  if (length(TFs) == 0) {
    message = paste0("No files with pattern ", par.l$TFBSPattern, " found in directory ", TFBS_dir, " as specified by the parameter \"dir_TFBS\"")
    checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
  }
  
  allTFs = sapply(strsplit(basename(TFs), "_"), "[[", 1)
  stopifnot(length(allTFs) > 0)
  
}


for (TFCur in allTFs) {
  
  TFCur = gsub(pattern = " ",replacement = "",TFCur)
  flog.info(paste0(" Checking TF ", TFCur, "..."))
  fileCur = paste0(TFBS_dir, "/", TFCur, par.l$TFBSPattern)
  if (!file.exists(fileCur)) {
    message = paste0("File ", fileCur, " does not exist even though the TF ", TFCur, " has been specified")
    checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
  }
  
  tableCur.df = read_tidyverse_wrapper(fileCur, type = "tsv", ncolExpected = 6, col_names = FALSE, col_types = "ciicnc")

  assertIntegerish(tableCur.df$X2, lower = 0)
  assertIntegerish(tableCur.df$X3, lower = 0)
  
  
}

flog.info(paste0("Done checking."))

# Write to dummy file
cat("DONE", file = par.l$output_flag)

.printExecutionTime(start.time)

flog.info("Session info: ", sessionInfo(), capture = TRUE)