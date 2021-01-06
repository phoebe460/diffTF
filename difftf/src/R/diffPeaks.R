start.time  <-  Sys.time()


#########################
# LIBRARY AND FUNCTIONS #
#########################

# Use the following line to load the Snakemake object to manually rerun this script (e.g., for debugging purposes)
# Replace {outputFolder} correspondingly.
# snakemake = readRDS("{outputFolder}/LOGS_AND_BENCHMARKS/diffPeaks.R.rds")

library("checkmate")
assertClass(snakemake, "Snakemake")
assertDirectoryExists(snakemake@config$par_general$dir_scripts)
source(paste0(snakemake@config$par_general$dir_scripts, "/functions.R"))

########################################################################
# SAVE SNAKEMAKE S4 OBJECT THAT IS PASSED ALONG FOR DEBUGGING PURPOSES #
########################################################################


createDebugFile(snakemake)

initFunctionsScript(packagesReq = NULL, minRVersion = "3.1.0", warningsLevel = 1, disableScientificNotation = TRUE)
checkAndLoadPackages(c("tidyverse", "futile.logger", "DESeq2", "csaw", "checkmate", "limma", "tools", "matrixStats"), verbose = FALSE)



###################
#### PARAMETERS ###
###################

par.l = list()

par.l$pseudocountAddition = 1
par.l$verbose = TRUE
par.l$log_minlevel = "INFO"

#####################
# VERIFY PARAMETERS #
#####################

checkAndLogWarningsAndErrors(snakemake, checkmate::checkClass(snakemake, "Snakemake"))

## INPUT ##
checkAndLogWarningsAndErrors(snakemake@input, checkmate::checkList(snakemake@input, min.len = 1))
checkAndLogWarningsAndErrors(snakemake@input, checkmate::checkSubset(names(snakemake@input), c("", "sampleData", "BAMPeakoverlaps")))

par.l$file_input_sampleData = snakemake@input$sampleData
checkAndLogWarningsAndErrors(par.l$file_input_sampleData, checkmate::checkFileExists(par.l$file_input_sampleData, access = "r"))

  
par.l$file_input_peakOverlaps = snakemake@input$BAMPeakoverlaps

for (fileCur in par.l$files_input_TF_summary) {
  checkAndLogWarningsAndErrors(fileCur, checkmate::checkFileExists(fileCur, access = "r"))
}

## OUTPUT ##
checkAndLogWarningsAndErrors(snakemake@output, checkmate::checkList(snakemake@output, min.len = 1))
checkAndLogWarningsAndErrors(names(snakemake@output), checkmate::checkSubset(c("", "sampleDataR", "peakFile", "peaks_tsv", "condComp", "normFacs", "normCounts", "plots", "DESeqObj"), names(snakemake@output)))

par.l$file_output_metadata     = snakemake@output$sampleDataR
par.l$file_output_peaks        = snakemake@output$peakFile
par.l$file_output_peaksTSV     = snakemake@output$peaks_tsv
# par.l$file_output_peaksPermTSV = snakemake@output$peaksPerm_tsv
par.l$file_output_condComp     = snakemake@output$condComp  
par.l$file_output_normFacs     = snakemake@output$normFacs
par.l$file_output_normCounts   = snakemake@output$normCounts
par.l$file_output_plots        = snakemake@output$plots
par.l$file_output_DESeqObj     = snakemake@output$DESeqObj


## CONFIG ##
checkAndLogWarningsAndErrors(snakemake@config,checkmate::checkList(snakemake@config, min.len = 1))

par.l$designFormula = snakemake@config$par_general$designContrast
checkAndLogWarningsAndErrors(par.l$designFormula, checkCharacter(par.l$designFormula, len = 1, min.chars = 3))

par.l$designFormulaVariableTypes = snakemake@config$par_general$designVariableTypes
checkAndLogWarningsAndErrors(par.l$designFormulaVariableTypes, checkCharacter(par.l$designFormulaVariableTypes, len = 1, min.chars = 3))

par.l$nPermutations = snakemake@config$par_general$nPermutations
checkAndLogWarningsAndErrors(par.l$nPermutations, checkIntegerish(par.l$nPermutations, lower = 0))

par.l$conditionComparison  = snakemake@config$par_general$conditionComparison
checkAndLogWarningsAndErrors(par.l$conditionComparison, checkCharacter(par.l$conditionComparison, len = 1))

par.l$debugMode = setDebugMode(snakemake@config$par_general$debugMode)

## PARAMS ##
checkAndLogWarningsAndErrors(snakemake@params,checkmate::checkList(snakemake@params, min.len = 1))
checkAndLogWarningsAndErrors(names(snakemake@params), checkSubset(names(snakemake@params), c("", "doCyclicLoess", "debugFile")))

par.l$doCyclicLoess = as.logical(snakemake@params$doCyclicLoess)
checkAndLogWarningsAndErrors(par.l$doCyclicLoess, checkFlag(par.l$doCyclicLoess))

## LOG ##
checkAndLogWarningsAndErrors(snakemake@log,checkmate::checkList(snakemake@log, min.len = 1))
par.l$file_log = snakemake@log[[1]]


allDirs = c(dirname(par.l$file_output_metadata), 
            dirname(par.l$file_output_peaks), 
            dirname(par.l$file_output_normFacs), 
            dirname(par.l$file_output_peaksTSV),
            dirname(par.l$file_log)
            )

testExistanceAndCreateDirectoriesRecursively(allDirs)


######################
# FINAL PREPARATIONS #
######################
startLogger(par.l$file_log, par.l$log_minlevel,  removeOldLog = TRUE)
printParametersLog(par.l)



if (par.l$debugMode) {
    
    flog.info(paste0("Debug mode active. Reading it all files that this step requires and save it to ", snakemake@params$debugFile))
    
    
    # Read all files already here and then save the session so as much as possible from the script can be executed without file dependencies
    sampleData.df = read_tidyverse_wrapper(par.l$file_input_sampleData, type = "tsv", col_names = TRUE, col_types = cols())
    coverageAll.df = read_tidyverse_wrapper(par.l$file_input_peakOverlaps, type = "tsv", col_names = TRUE, comment = "#", col_types = cols())
    
    save(list = ls(), file = snakemake@params$debugFile)
    flog.info(paste0("File ", snakemake@params$debugFile, " has been saved. You may use it for trouble-shooting and debugging, see the Documentation for more details."))
    
}


#################
# READ METADATA #
#################

sampleData.df = read_tidyverse_wrapper(par.l$file_input_sampleData, type = "tsv", col_names = TRUE, col_types = cols())

checkAndLogWarningsAndErrors(colnames(sampleData.df), checkSubset(c("bamReads"), colnames(sampleData.df)))

designComponents.l = checkDesignIntegrity(snakemake, par.l, sampleData.df)

components3types   = designComponents.l$types
variableToPermute  = designComponents.l$variableToPermute

# Which of the two modes should be done, pairwise or quantitative?
comparisonMode = "quantitative"
if (components3types["conditionSummary"] == "logical" | components3types["conditionSummary"] == "factor") {
    comparisonMode = "pairwise"
}

# Read and modify samples metadata
sampleData.df = mutate(sampleData.df, name = file_path_sans_ext(basename(sampleData.df$bamReads)),
                       SampleID = as.character(SampleID))
  

# Check and change column types as specified in the design formula
for (colnameCur in names(components3types)) {
  
  coltype = components3types[colnameCur]
  if (coltype == "factor") {
    sampleData.df[,colnameCur] = as.factor(unlist(sampleData.df[,colnameCur]))
  } else if (coltype == "numeric") {
    sampleData.df[,colnameCur] = as.numeric(unlist(sampleData.df[,colnameCur]))
  } else if (coltype == "integer") {
    sampleData.df[,colnameCur] = as.integer(unlist(sampleData.df[,colnameCur]))
  } else if (coltype == "logical") {
    sampleData.df[,colnameCur] = as.logical(unlist(sampleData.df[,colnameCur]))
  } else {
    message = paste0("Unknown type: ", colnameCur)
    checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
  }
  
}

if (comparisonMode == "pairwise") {
    
    # Change the conditionSummary specifically and enforce the direction as specified in the config file
    sampleData.df$conditionSummary = factor(sampleData.df$conditionSummary, levels = strsplit(par.l$conditionComparison, ",")[[1]])
} 


##############################
# ITERATE THROUGH PEAK FILES #
##############################

coverageAll.df = read_tidyverse_wrapper(par.l$file_input_peakOverlaps, type = "tsv",
                                        col_names = TRUE, comment = "#", col_types = cols())

## transform as matrix data frame with counts
coverageAll.m = as.matrix(dplyr::select(coverageAll.df, -one_of("Geneid", "Chr", "Start", "End", "Strand", "Length")))

# Take the basenames of the files, which have not been modified in the sorted versions of the BAMs as they are just located in different folders
sampleIDs = sampleData.df$SampleID[which(basename(sampleData.df$bamReads) %in% basename(colnames(coverageAll.m)))]

if (length(unique(sampleIDs)) != nrow(sampleData.df)) {
    message = paste0("Colnames mismatch. Make sure that each sampleID is unique in the sample summary table.")
    checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
} 

if (length(sampleIDs) != ncol(coverageAll.m)) {
  message = paste0("Mismatch between number of sample IDs (", length(sampleIDs), ") and number of columns in coverage file (", ncol(coverageAll.m), "). It appears that the number of samples been changed after running the pipeline the first time. Rerun the full pipeline from scratch.")
  checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
} 

colnames(coverageAll.m) = sampleIDs

rownames(coverageAll.m) = coverageAll.df$Geneid


peaks.df = dplyr::select(coverageAll.df, one_of(c("Chr", "Start", "End", "Geneid")))
colnames(peaks.df) = c("chr", "PSS", "PES", "annotation")
# Filter and retain only unique identifiers
peaks.filtered.df = distinct(peaks.df, annotation, .keep_all = TRUE)
nRowsFiltered = nrow(peaks.df) - nrow(peaks.filtered.df)
if (par.l$verbose & nRowsFiltered  > 0) flog.info(paste0("Filtered ", nRowsFiltered, " non-unique positions out of ", nrow(peaks.df), " from peaks table."))
peaks.df = peaks.filtered.df

# Save the last, they are all identical anyway except for the count column
saveRDS(peaks.filtered.df, file = par.l$file_output_peaks)


#############################
# SAMPLE LABEL PERMUTATIONS #
#############################

sampleData.l = list()
nSamples = nrow(sampleData.df)
sampleData.l[["permutation0"]] = sampleData.df

conditionCounter = table(sampleData.df[,variableToPermute])



if (comparisonMode == "pairwise") {
    
    # Record the frequency of the conditions to determine how many permutations are possible
    nSamplesRareCondition     = min(conditionCounter)
    nSamplesFrequentCondition = max(conditionCounter)
    nameRareCondition         = names(conditionCounter)[conditionCounter == min(conditionCounter)][1]
    nameFrequentCondition     = names(conditionCounter)[which(names(conditionCounter) != nameRareCondition)]
    nPermutationsTotal        = choose(nSamples, nSamplesFrequentCondition) # same as choose(nSamples, nSamplesRareCondition)
  
} else {
    
    # The total number of permutations is calculated using "Permutations of multisets"
    # nPermutationsTotal = factorial(nSamples)/product(factorial(conditionCounter))
    # Alternatively, calculate as choose(nSamples, groupSize1)*choose(nSamples-groupSize1, groupSize2)* ... * choose(nSamples - groupSize1 -groupSize2 - ... - groupSize(n-2), groupSize(n-1))
    # Generate a vector of the first components of the choose argument
    productVec = nSamples - sapply(1:(length(conditionCounter)-1), function(x) {sum(conditionCounter[1:x])})
    nPermutationsTotal = matrixStats::product(choose(c(nSamples, productVec),conditionCounter[-length(conditionCounter)])) # might be better because it does not calculate nSamples in the first place
    
}

if (nPermutationsTotal < par.l$nPermutations) {
    
    message = paste0("The total number of possible permutations is only ", nPermutationsTotal, ", but more have been requested. The value for the parameter nPermutations will be adjusted.")
    checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
    par.l$nPermutations = nPermutationsTotal
}
    
# Permute samples beforehand here so that each call to a permutation is unique
permutationsList.l = list()
nPermutationsDone = 0
failsafeCounter   = 0
while (nPermutationsDone < par.l$nPermutations) {
    
    # TODO: Do not include the original, non-shuffled variant here
    sampleCur = sample.int(nSamples)
    
    # Check whether this is a "new" permutation 
    if (comparisonMode == "pairwise") {
        
        samplesRareCondShuffled  = sampleData.df$SampleID[which(sampleData.df$conditionSummary[sampleCur] == nameRareCondition)]
        indexNameCur = paste0(sort(samplesRareCondShuffled), collapse = ",")
        
    } else {
        
        # Slighlty more complicated here because we may have more than two groups, do it per group then
        
        indexNameCur = ""
        for (groupCur in 1:length(conditionCounter)) {
            samplesRareCondShuffled  = sampleData.df$SampleID[which(sampleData.df$conditionSummary[sampleCur] == as.numeric(names(conditionCounter)[groupCur]))]
            indexNameCurGroup = paste0(sort(samplesRareCondShuffled), collapse = ",")
            indexNameCur = paste0(indexNameCur, indexNameCurGroup, "+")
        }
        
    }
    
    
    # Check if this permutation has already been used. If yes, produce a different one
    
    if (!indexNameCur %in% names(permutationsList.l)) {
        failsafeCounter   = 0
        permutationsList.l[[indexNameCur]] = sampleCur
        nPermutationsDone = nPermutationsDone + 1
    } else {

        failsafeCounter   =  failsafeCounter + 1
        if (failsafeCounter > 5000) {
            message = "Could not generate more permutations. This looks like a bug."
            checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
        } 
    }
    
} 


##############################################
# RUN DESEQ TO OBTAIN NORMALIZED COUNTS ONLY #
##############################################

designFormula = convertToFormula(par.l$designFormula, colnames(sampleData.df))


# Enforce the correct order of sampleData and coverageAll so that the rownames match
stopifnot(sampleData.df$SampleID %in% colnames(coverageAll.m))
sampleData.temp.df = as.data.frame(sampleData.df)
rownames(sampleData.temp.df) = sampleData.temp.df$SampleID

cds.peaks <- DESeqDataSetFromMatrix(countData = coverageAll.m[, sampleData.temp.df$SampleID],
                                    colData = sampleData.temp.df,
                                    design = designFormula)

# Recent versions of DeSeq seem to do this automatically, whereas older versions don't, so enforce it here
if (!identical(colnames(cds.peaks), colnames(coverageAll.m))) {
    colnames(cds.peaks) = colnames(coverageAll.m)
}
if (!identical(rownames(cds.peaks), rownames(coverageAll.m))) {
    rownames(cds.peaks) = rownames(coverageAll.m)
}

# Do a regular size factor normalization
if (!par.l$doCyclicLoess) {
    
    cds.peaks <- estimateSizeFactors(cds.peaks)
    
    normFacs = sizeFactors(cds.peaks)
    
    # TODO: NormFacts identical between different permutations? yes or?
    
} else {
    
    # Perform a cyclic loess normalization
    # We use a slighlty more complicated setup to derive size factors for library normalization
    # Instead of just determining the size factors in DeSeq2 via cirtual samples, we use 
    # a normalization from the csaw package (see https://www.rdocumentation.org/packages/csaw/versions/1.6.1/topics/normOffsets)
    # and apply a non-linear normalization. 
    # For each sample, a lowess curve is fitted to the log-counts against the log-average count. 
    # The fitted value for each bin pair is used as the generalized linear model offset for that sample. 
    # The use of the average count provides more stability than the average log-count when low counts are present for differentially bound regions.
    
    
    # since counts returns,by default, non-normalized counts, the following code should be fine and there is no need to
    # also run estimateSizeFactors beforehand
    
    if (packageVersion("csaw") <= "1.14.1") {
        normFacs = exp(normOffsets(DESeq2::counts(cds.peaks), lib.sizes = colSums(DESeq2::counts(cds.peaks)), type = "loess"))
    } else {
        object = SummarizedExperiment(list(counts=DESeq2::counts(cds.peaks)))
        object$totals = colSums(DESeq2::counts(cds.peaks))
        normFacs  = exp(normOffsets(object, se.out = FALSE))
    }

    
    # sanity check: is the geometric mean across samples equal to one?
    #library("psych")
    #all.equal(geometric.mean(t(normFacs)), rep(1, dim(cds.peaks)[1]))
    
    rownames(normFacs) = rownames(coverageAll.m)
    colnames(normFacs) = colnames(coverageAll.m)
    
    # We now provide gene-specific normalization factors for each sample as a matrix, which will preempt sizeFactors
    normalizationFactors(cds.peaks) <- normFacs
    
}


# Filter peaks with zero counts
cds.peaks.filt   = cds.peaks[rowMeans(DESeq2::counts(cds.peaks)) > 0, ]


# DESeq log2fc are not used at all afterwards, as we currently only take the normalization factors to normalize the TFBS subsequently


if (comparisonMode == "pairwise") {
    
    # The levels have to be reversed because the first element is the one appearing at the right of the plot, with positive values as. compared to the reference
    comparisonDESeq = rev(levels(sampleData.df$conditionSummary))
    
} else {
    
    comparisonDESeq = c("positive change", "negative change")
}


##############
# GET LOG2FC #
##############

countsNorm        = DESeq2::counts(cds.peaks.filt, norm = TRUE)
countsNorm.df     = as.data.frame(countsNorm) %>%
  dplyr::mutate(peakID = rownames(cds.peaks.filt))  %>%
  dplyr::select(one_of("peakID", colnames(countsNorm)))


# Deseq analysis
cds.peaks.filt = tryCatch( {
    DESeq(cds.peaks.filt, fitType = 'local', quiet = TRUE)

}, error = function(e) {
    message = "Warning: Could not run DESeq with local fitting, retry with default fitting type..."
    checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
    DESeq(cds.peaks.filt, quiet = TRUE)
}
)
#Enforce the correct order of the comparison
if (comparisonMode == "pairwise") { 
  
  contrast = c(variableToPermute, comparisonDESeq[1], comparisonDESeq[2])

} else {
  # Same as without specifying contrast at all
  contrast = list(variableToPermute)
     
}
cds.peaks.df <- as.data.frame(DESeq2::results(cds.peaks.filt, contrast = contrast)) 


final.peaks.df = tibble( 
    "permutation" = 0,
    "peakID"    = rownames(cds.peaks.df), 
    "DESeq_baseMean" = cds.peaks.df$baseMean,
    "l2FC"     = cds.peaks.df$log2FoldChange,
    "DESeq_ldcSE"    = cds.peaks.df$lfcSE,
    "DESeq_stat"     = cds.peaks.df$stat,
    "pval"     =  cds.peaks.df$pvalue, 
    "pval_adj" =  cds.peaks.df$padj
)

plotDiagnosticPlots(cds.peaks.filt, comparisonDESeq, file = par.l$file_output_plots, contrast = contrast, maxPairwiseComparisons = 20)


saveRDS(cds.peaks.filt, file = par.l$file_output_DESeqObj)

####################
# RUN PERMUTATIONS #
####################

if (par.l$nPermutations > 0) {
  #Rename so it is easier to address in the following code
  listNames = paste0("permutation", seq_len(par.l$nPermutations))
  names(permutationsList.l) = listNames
  
  # final.peaks.perm.df = tribble(~permutation, ~peakID, ~l2FC)
  sampleDataOrig.df = sampleData.df
  
  for (permutationCur in names(permutationsList.l)) {
    
    # flog.info(paste0("Running for permutation ", permutationCur))
    
    sampleData.df = sampleDataOrig.df
    sampleData.df[,variableToPermute] = unlist(sampleData.df[,variableToPermute]) [permutationsList.l[[permutationCur]]]
    
    sampleData.l[[permutationCur]] = sampleData.df

    # We don't need permuted peaks l2fc actually so permutation-specific l2fc can be skipped
  
  }
}



################
# WRITE OUTPUT #
################

saveRDS(sampleData.l, par.l$file_output_metadata)

# Do it separately for each column because different rounding schemes might be needed
final.peaks.df = final.peaks.df %>% 
    mutate(permutation = as.integer(permutation)) %>%
    mutate_at(c("DESeq_baseMean", "l2FC", "DESeq_ldcSE", "DESeq_stat"), signif, 3) %>%
    mutate_at(c("pval", "pval_adj"), formatC, format = "g", digits = 3)


write_tsv(final.peaks.df, path = par.l$file_output_peaksTSV)

#final.peaks.perm.df = mutate_if(final.peaks.perm.df, is.numeric, as.character)
#write_tsv(final.peaks.perm.df, path = par.l$file_output_peaksPermTSV)

saveRDS(comparisonDESeq, file = par.l$file_output_condComp)

saveRDS(normFacs, par.l$file_output_normFacs)

# Deactivated, because column names starting with numbers will cause problems and crash. Since the names come from the sample file, this cannot be excluded.
#countsNorm.df = mutate_if(countsNorm.df, is.numeric, as.character)
write_tsv(countsNorm.df, path = par.l$file_output_normCounts)



.printExecutionTime(start.time)

flog.info("Session info: ", sessionInfo(), capture = TRUE)
