start.time  <-  Sys.time()

#########################
# LIBRARY AND FUNCTIONS #
#########################
library("checkmate")
assertClass(snakemake, "Snakemake")
assertDirectoryExists(snakemake@config$par_general$dir_scripts)
source(paste0(snakemake@config$par_general$dir_scripts, "/functions.R"))

########################################################################
# SAVE SNAKEMAKE S4 OBJECT THAT IS PASSED ALONG FOR DEBUGGING PURPOSES #
########################################################################

# Use the following line to load the Snakemake object to manually rerun this script (e.g., for debugging purposes)
# Replace {outputFolder} correspondingly.
# snakemake = readRDS("{outputFolder}/LOGS_AND_BENCHMARKS/summaryFinal.R.rds")
# Note that one currently cannot overwrite valeus within the snakemake object; instead, assign them to a variable as done below and change if necessary
createDebugFile(snakemake)

initFunctionsScript(packagesReq = NULL, minRVersion = "3.1.0", warningsLevel = 1, disableScientificNotation = TRUE)
checkAndLoadPackages(c("tidyverse", "futile.logger", "ggrepel", "checkmate", "tools", "grDevices", "locfdr"), verbose = FALSE)


###################
#### PARAMETERS ###
###################

# Currently hard-coded parameters
par.l = list()

# 1. Misc
par.l$verbose = TRUE
par.l$log_minlevel = "INFO"

# 2. Statistical thresholds and values
par.l$significanceThresholds  = c(0.001, 0.01, 0.05,0.1,0.2) # p-value thresholds
par.l$classes_CohensD = c("small", "medium", "large", "very large")
par.l$thresholds_CohensD = c(0.1, 0.5, 0.8)

# 3. RNA-Seq specific
par.l$corMethod = "pearson" # Expression-peak count correlation method. As we quantile normalize now, should be pearson
par.l$regressionMethod = "glm" # for correlating RNA-Seq classification with TF activity for function plot_classCorrelations
par.l$normMethodRNA = "quantile"
par.l$filter_minCountsPerCondition = 5 # For filtering RNA-seq genes, see Documentation
par.l$filter_minMedianAll = 0 # For filtering RNA-seq genes, see Documentation 
par.l$thresholds_pvalue_Wilcoxon = 0.05
# Stringencies for AR classification
par.l$allClassificationThresholds = c(0.1, 0.05, 0.01, 0.001)

# 4. Volcano plot settings
par.l$maxTFsToLabel = 150 # Maximum Tfs to label in the Volcano plot
par.l$volcanoPlot_minDimensions  = 12
par.l$minPointSize = 0.3
par.l$plot_grayColor = "grey50"
par.l$colorCategories = c("activator" = "#4daf4a", "undetermined" = "black", "repressor" = "#e41a1c", "not-expressed" = "Snow3") # diverging, modified
par.l$colorConditions = c("#ef8a62", "#67a9cf") # Colors of the two conditions for the background

#####################
# VERIFY PARAMETERS #
#####################

assertClass(snakemake, "Snakemake")

## INPUT ##
assertList(snakemake@input, min.len = 1)
assertSubset(c("", "allPermutationResults", "condComp", "normCounts"), names(snakemake@input))


par.l$files_input_permResults  = snakemake@input$allPermutationResults
for (fileCur in par.l$files_input_permResults) {
    assertFileExists(fileCur, access = "r")
}

par.l$file_input_condCompDeSeq = snakemake@input$condComp
assertFileExists(par.l$file_input_condCompDeSeq, access = "r")

par.l$file_input_countsNorm = snakemake@input$normCounts
assertFileExists(par.l$file_input_countsNorm, access = "r")

par.l$file_input_metadata = snakemake@input$sampleDataR
assertFileExists(par.l$file_input_metadata, access = "r")

## OUTPUT ##
assertList(snakemake@output, min.len = 1)
assertSubset(c("", "summary", "volcanoPlot", "diagnosticPlots", "plotsRDS"), names(snakemake@output))

par.l$file_output_summary  = snakemake@output$summary
par.l$file_plotVolcano     = snakemake@output$volcanoPlot
par.l$files_plotDiagnostic = snakemake@output$diagnosticPlots
par.l$file_output_plots    = snakemake@output$plotsRDS

## CONFIG ##
assertList(snakemake@config, min.len = 1)

par.l$plotRNASeqClassification = as.logical(snakemake@config$par_general$RNASeqIntegration)
assertFlag(par.l$plotRNASeqClassification)

par.l$nPermutations = snakemake@config$par_general$nPermutations
assertIntegerish(par.l$nPermutations, lower = 0)

par.l$outdir = snakemake@config$par_general$outdir

par.l$debugMode = setDebugMode(snakemake@config$par_general$debugMode)

if (par.l$plotRNASeqClassification) {
    
    par.l$file_input_HOCOMOCO_mapping    = snakemake@config$additionalInputFiles$HOCOMOCO_mapping
    par.l$file_input_geneCountsPerSample = snakemake@config$additionalInputFiles$RNASeqCounts
    assertFileExists(par.l$file_input_HOCOMOCO_mapping, access = "r")
    assertFileExists(par.l$file_input_geneCountsPerSample, access = "r")
}

## LOG ##
assertList(snakemake@log, min.len = 1)
par.l$file_log = snakemake@log[[1]]

allDirs = c(dirname(par.l$file_output_summary), dirname(par.l$files_plotDiagnostic),dirname(par.l$file_log))
testExistanceAndCreateDirectoriesRecursively(allDirs)

assertCharacter(par.l$colorCategories, len = 4)
assertSubset(names(par.l$colorCategories), c("activator", "undetermined", "repressor", "not-expressed"))
assertCharacter(par.l$colorConditions, len = 2)

######################
# FINAL PREPARATIONS #
######################
startLogger(par.l$file_log, par.l$log_minlevel,  removeOldLog = TRUE)
printParametersLog(par.l)

if (par.l$debugMode) {
    
    flog.info(paste0("Debug mode active. Reading it all files that this step requires and save it to ", snakemake@params$debugFile))
    
    # Read all files already here and then save the session so as much as possible from the script can be executed without file dependencies
    conditionComparison = readRDS(par.l$file_input_condCompDeSeq)
    
    results.l = list()
    for (fileCur in par.l$files_input_permResults) {
        results.l[[fileCur]] =  read_tidyverse_wrapper(fileCur, type = "tsv", col_names = TRUE, col_types = cols())
    }
    output.global.TFs.orig = do.call(rbind.data.frame, results.l)

    # Remove rows with NA
    TF_NA = which(is.na(output.global.TFs.orig$weighted_meanDifference))
    if (length(TF_NA) > 0) {
        output.global.TFs.orig = output.global.TFs.orig[-TF_NA,]
    }
    
    sampleData.l = readRDS(par.l$file_input_metadata)
    
    
    if (par.l$plotRNASeqClassification) {
        
        rootOutdir = snakemake@config$par_general$outdir
        assertCharacter(rootOutdir)
        
        comparisonType = snakemake@config$par_general$comparisonType
        assertCharacter(comparisonType)
        
        if (nchar(comparisonType) > 0) {
            comparisonType = paste0(comparisonType, ".")
        }
        
        extensionSize = as.integer(snakemake@config$par_general$regionExtension)
        assertIntegerish(extensionSize)

        countsRNA.all.df = read_tidyverse_wrapper(par.l$file_input_geneCountsPerSample, type = "tsv", col_names = TRUE)
        HOCOMOCO_mapping.df = readHOCOMOCOTable(par.l$file_input_HOCOMOCO_mapping)
        countsATAC.norm.df = read_tidyverse_wrapper(par.l$file_input_countsNorm, type = "tsv", col_types = cols())

        HOCOMOCO_mapping.df.overlap = filterHOCOMOCOTable(HOCOMOCO_mapping.df, unique(output.global.TFs.orig$TF))
        
        TF.peakMatrix.df = createBindingMatrixFromFiles(HOCOMOCO_mapping.df.overlap, countsATAC.norm.df, rootOutdir, extensionSize, comparisonType)
    }
    
    
    save(list = ls(), file = snakemake@params$debugFile)
    flog.info(paste0("File ", snakemake@params$debugFile, " has been saved. You may use it for trouble-shooting and debugging, see the Documentation for more details."))
}


################
# COLLECT DATA #
################

conditionComparison = readRDS(par.l$file_input_condCompDeSeq)
assertVector(conditionComparison, len = 2)

# Assemble the final table and collect permutation information from all TFs
output.global.TFs.orig = NULL

nTF = length(par.l$files_input_permResults)
for (fileCur in par.l$files_input_permResults) {
    
    resultsCur.df =  read_tidyverse_wrapper(fileCur, type = "tsv", col_names = TRUE, col_types = cols())

    assertIntegerish(nrow(resultsCur.df), lower = 1, upper = par.l$nPermutations + 1)
    
    if (is.null(output.global.TFs.orig)) {
        output.global.TFs.orig = resultsCur.df
    } else {
        output.global.TFs.orig = rbind(output.global.TFs.orig, resultsCur.df)
    }
    
}

# Convert columns to numeric if they are not already
output.global.TFs.orig = mutate(output.global.TFs.orig,
                                weighted_meanDifference = as.numeric(weighted_meanDifference),
                                variance                = as.numeric(variance),
                                weighted_CD             = as.numeric(weighted_CD),
                                weighted_Tstat          = as.numeric(weighted_Tstat))

# Remove rows with NA
TF_NA = which(is.na(output.global.TFs.orig$weighted_meanDifference))
permNACount = table(output.global.TFs.orig$permutation[TF_NA])

if (length(TF_NA) > 0) {
    
    permWithNA = names(permNACount)
    nPermOnlyNA = length(which(permNACount == nTF))
    if (nPermOnlyNA > 0) {
        message = paste0("Data from ",  nPermOnlyNA, " permutations had to be removed due to NA values in weighted_meanDifference (insufficient data in previous steps). This frequently happens for small datasets and only affects the calculation of p-values. ")
        checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
    }
    
    if ("0" %in% permWithNA) {
        message = paste0("Data from ",  permNACount["0"], " TF had to be removed due to NA values in weighted_meanDifference (insufficient data in previous steps).")
        checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
    }
    
    flog.info(paste0("Removing the following TFs from the output table due to missing data: ", paste0(output.global.TFs.orig$TF[TF_NA], collapse = ",")))
    output.global.TFs.orig = output.global.TFs.orig[-TF_NA,]
    
    
}


# Make sure to only use as many permutations as actually could be done, independent of what the user specified before
par.l$nPermutations = max(output.global.TFs.orig$permutation)

########################################################
# FILTER BY PERMUTATIONS AND COMPARE, DIAGNOSTIC PLOTS #
########################################################
# Compare the distributions from the real and random permutations
diagPlots.l = list()

pdf(par.l$files_plotDiagnostic[1])
output.global.TFs.orig$pvalue = NA
if (par.l$nPermutations > 0) {
    
    dataReal.df = dplyr::filter(output.global.TFs.orig, permutation == 0)
    dataPerm.df = dplyr::filter(output.global.TFs.orig, permutation > 0)
    
    
    # ( (#perm >TH)/#perm ) / ( (#perm >TH)/#perm + (#real > TH)/n_real )
    
    xrange = range(output.global.TFs.orig$weighted_meanDifference, na.rm = TRUE)
    
    maxY = max(density(dataPerm.df$weighted_meanDifference, na.rm = TRUE)$y, density(dataReal.df$weighted_meanDifference, na.rm = TRUE)$y)
    
    plot(density(dataPerm.df$weighted_meanDifference, na.rm = TRUE), col = "black", main = "Weighted mean difference values (black = permuted)", 
         xlim = xrange, ylim = c(0,maxY))
    lines(density(dataReal.df$weighted_meanDifference, na.rm = TRUE), col = "red")
    
    for (TFCur in unique(output.global.TFs.orig$TF)) {
        
        dataCur.df = dplyr::filter(output.global.TFs.orig, TF == TFCur)
        dataReal.df = dplyr::filter(dataCur.df, permutation == 0)
        dataPerm.df = dplyr::filter(dataCur.df, permutation > 0)
        
        
        if (nrow(dataPerm.df) == 0) {
            message = paste0("For TF : ", TFCur, ", no permutation data could be found. This is not supposed to happen. Rerun the binning step for this TF.")
            checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
            next
        }
        
        rangeX = range(dataCur.df$weighted_meanDifference)
        rowCur = which(output.global.TFs.orig$TF == TFCur)
        
        g = ggplot(dataPerm.df, aes(weighted_meanDifference)) + geom_density() + 
            geom_vline(xintercept = dataReal.df$weighted_meanDifference[1], color = "red") + ggtitle(TFCur) + xlim(c(rangeX * c(1.5,1.5))) 
        # + scale_x_continuous(limits = c(min(dataCur.df$weighted_meanDifference) - 0.5, max(dataCur.df$weighted_meanDifference) + 0.5))
        
        diagPlots.l[[TFCur]] = g
        plot(g)
        
        nPermThreshold = length(which(abs(dataPerm.df$weighted_meanDifference) > abs(dataReal.df$weighted_meanDifference[1])))
        
        # We add a 1 to the numerator and denominator to account for misestimation of the p-value 
        # (for more details see Phipson and Smyth, Permutation P-values should never be zero).
        # See here: https://www.ncbi.nlm.nih.gov/pubmed/21044043
        pvalueCur = (nPermThreshold + 1) / (par.l$nPermutations + 1)
        output.global.TFs.orig$pvalue[rowCur] = pvalueCur 
        
    }
    
    # TODO: Diagnostic plot for pvalue
    plot(ggplot(output.global.TFs.orig, aes(pvalue)) + geom_density() + ggtitle("Local fdr density across all TF"))
    
} else {
    
    populationMean = 0
    zScore = (output.global.TFs.orig$weighted_Tstat - populationMean) / sqrt(output.global.TFs.orig$variance)
    
    # 2-sided test
    output.global.TFs.orig$pvalue   = 2*pnorm(-abs(zScore))
    
    # Handle extreme cases with p-values that are practically 0 and would cause subsequent issues
    index0 = which(output.global.TFs.orig$pvalue < .Machine$double.xmin)
    if (length(index0) > 0) {
        output.global.TFs.orig$pvalue[index0] = .Machine$double.xmin
    }
}

output.global.TFs.permutations = dplyr::filter(output.global.TFs.orig, permutation > 0)
output.global.TFs              = dplyr::filter(output.global.TFs.orig, permutation == 0)



output.global.TFs = mutate(output.global.TFs, 
                           Cohend_factor = case_when( weighted_CD < par.l$thresholds_CohensD[1] ~ par.l$classes_CohensD[1],
                                                      weighted_CD < par.l$thresholds_CohensD[2] ~ par.l$classes_CohensD[2],
                                                      weighted_CD < par.l$thresholds_CohensD[3] ~ par.l$classes_CohensD[3],
                                                      TRUE ~ par.l$classes_CohensD[4]),
                           Cohend_factor = factor(Cohend_factor, levels = par.l$classes_CohensD, labels = seq_len(length(par.l$classes_CohensD))),
                           pvalueAdj     = p.adjust(pvalue, method = "BH"),
                           weighted_meanDifference = as.numeric(weighted_meanDifference))


colnamesToPlot = c("weighted_meanDifference", "weighted_CD", "TFBS", "weighted_Tstat", "variance", "pvalue", "pvalueAdj")

for (pValueCur in c(par.l$significanceThresholds , 1)) {
    
    filtered.df = dplyr::filter(output.global.TFs, pvalueAdj <= pValueCur)
    
    title = paste0("p-value: ", pValueCur, " (retaining ", nrow(filtered.df), " TF)")
    
    for (measureCur in colnamesToPlot) {
        
        if (all(!is.finite(unlist(filtered.df[,measureCur])))) {
            next
        }
        
        if (! measureCur %in% colnames(output.global.TFs)) {
            next
        }
        
        if (measureCur %in%  c("Cohend_factor")) {
            plot(ggplot(filtered.df, aes_string(measureCur))  + stat_count() + theme_bw() + ggtitle(title))
            
        } else {
            plot(ggplot(filtered.df, aes_string(measureCur))  + geom_histogram(bins = 50) + theme_bw() + ggtitle(title))
        }
        
    }
}

stats.df = group_by(output.global.TFs.orig, permutation) %>% summarise(max = max(weighted_meanDifference), min = min(weighted_meanDifference))
ggplot(stats.df, aes(min)) + geom_density()
ggplot(stats.df, aes(max)) + geom_density()

dev.off()


#################################
# mode of change quantification #
#################################

sampleData.l = readRDS(par.l$file_input_metadata)
sampleData.df = sampleData.l[["permutation0"]]
designComponents.l = checkDesignIntegrity(snakemake, par.l, sampleData.df)

components3types   = designComponents.l$types

# Which of the two modes should be done, pairwise or quantitative?
comparisonMode = "quantitative"
if (components3types["conditionSummary"] == "logical" | components3types["conditionSummary"] == "factor") {
    comparisonMode = "pairwise"
}


##########################
# INTEGRATE RNA-Seq DATA #
##########################
if (par.l$plotRNASeqClassification) {
    
    # Require some more packages here
    checkAndLoadPackages(c( "lsr", "DESeq2",  "matrixStats",  "pheatmap", "preprocessCore"), verbose = FALSE)
    
    
    classesList.l = list(c("activator","undetermined","repressor","not-expressed"),
                         c("activator","undetermined","repressor"),
                         c("activator","repressor")
    )
    
    extensionSize = as.integer(snakemake@config$par_general$regionExtension)
    assertIntegerish(extensionSize)
    
    rootOutdir = snakemake@config$par_general$outdir
    assertCharacter(rootOutdir)
    
    comparisonType = snakemake@config$par_general$comparisonType
    assertCharacter(comparisonType)
    
    if (nchar(comparisonType) > 0) {
        comparisonType = paste0(comparisonType, ".")
    }
    
    ########################
    # Process RNA-Seq data #
    ########################
    countsRNA.all.df = read_tidyverse_wrapper(par.l$file_input_geneCountsPerSample, type = "tsv", col_names = TRUE)
 
    # The design formula for RNA-Seq is different from the one we used before for ATAC-Seq
    # Either take the one that the user provided or, if he did not, use a general one with only the condition
    # The formula is only needed for the function compute_classCorrelations below
    par.l$designFormulaRNA = snakemake@config$par_general$designContrastRNA
    if (is.null(par.l$designFormulaRNA)) {
        par.l$designFormulaRNA = "~conditionSummary"
        flog.warn(paste0("Could not find the parameter designContrastRNA in the configuration file. The default of \"~conditionSummary\" will be taken as formula. If you know about confounding variables, rerun this step and add the parameter (see the Documentation for details)"))
    }
    designFormulaRNA = convertToFormula(par.l$designFormulaRNA, colnames(sampleData.df))
    
    countsRNA.all.mod.df = as.data.frame(countsRNA.all.df[,sampleData.df$SampleID])
    rownames(countsRNA.all.mod.df) = gsub("\\..+", "", countsRNA.all.df$ENSEMBL, perl = TRUE)

    dd <- DESeqDataSetFromMatrix(countData = countsRNA.all.mod.df,
                                 colData = sampleData.df,
                                 design = designFormulaRNA)
    
    dd = estimateSizeFactors(dd)
    
    ######################################
    # Filtering of lowly expressed genes #
    ######################################
    dd.filt = filterLowlyExpressedGenes(dd, comparisonMode, par.l$filter_minCountsPerCondition, par.l$filter_minMedianAll)
    
    # Raw counts, used for other types of normalization thereafter
    dd.filt.rawCounts =  DESeq2::counts(dd.filt, normalized=FALSE)

    countsRNA.norm.df  = normalizeCounts(dd.filt.rawCounts, method = par.l$normMethodRNA, idColumn = NULL)
    
    # Loading TF gene translation table
    HOCOMOCO_mapping.df = readHOCOMOCOTable(par.l$file_input_HOCOMOCO_mapping)
    
    ####################
    # READ ATAC counts #
    ####################
    # Already normalized, coming from previous steps
    countsATAC.norm.df = read_tidyverse_wrapper(par.l$file_input_countsNorm, type = "tsv", col_types = cols())

    ##########################
    # Intersect RNA and ATAC #
    ########################## 
    
    # Subset data to retain only samples that appear in both RNA and ATAC
    data.l = intersectData(countsRNA.norm.df, countsATAC.norm.df)
    
    countsRNA.norm.df  = data.l[["RNA"]]
    countsATAC.norm.df = data.l[["ATAC"]] 
  
    # Filter genes that might not be in the output.global.TFs$TF list 
    HOCOMOCO_mapping.df.overlap = filterHOCOMOCOTable(HOCOMOCO_mapping.df, output.global.TFs$TF)
   
    TF.peakMatrix.df = createBindingMatrixFromFiles(HOCOMOCO_mapping.df.overlap, countsATAC.norm.df, rootOutdir, extensionSize, comparisonType)
    
    res.l = filterPeaksByRowMeans(countsATAC.norm.df, TF.peakMatrix.df, minMean = 1)
    TF.peakMatrix.df        = res.l[["bindingMatrix"]]  
    countsATAC.norm.filt.df = res.l[["peakCounts"]]  
    

    sort.cor.m = correlateATAC_RNA(countsRNA.norm.df, countsATAC.norm.filt.df, HOCOMOCO_mapping.df.overlap, corMethod = "pearson")
    
    res.l = computeForegroundAndBackgroundMatrices(TF.peakMatrix.df, sort.cor.m)
    median.cor.tfs       = res.l[["median_foreground"]]
    median.cor.tfs.non   = res.l[["median_background"]]
    t.cor.sel.matrix     = res.l[["foreground"]]
    t.cor.sel.matrix.non = res.l[["background"]]
    
    # 3. Final classification: Calculate thresholds by calculating the quantiles of the background anhd compare the real values to the background
    act.rep.thres.l = calculate_classificationThresholds(t.cor.sel.matrix.non, par.l)
    
    output.global.TFs = finalizeClassificationAndAppend(output.global.TFs, median.cor.tfs, act.rep.thres.l, par.l, t.cor.sel.matrix, t.cor.sel.matrix.non, significanceThreshold_Wilcoxon = par.l$thresholds_pvalue_Wilcoxon)
    
    ####################
    ####################
    # DIAGNOSTIC PLOTS #
    ####################
    ####################
    HOCOMOCO_mapping.df.overlap.exp = dplyr::filter(HOCOMOCO_mapping.df.overlap, ENSEMBL %in% countsRNA.norm.df$ENSEMBL)
    
    pdf(file = par.l$files_plotDiagnostic[2], width = 4, height = 8)
    plot_AR_thresholds(median.cor.tfs, median.cor.tfs.non, par.l, act.rep.thres.l, file = NULL)
    plot_heatmapAR(TF.peakMatrix.df, HOCOMOCO_mapping.df.overlap.exp, sort.cor.m, par.l, median.cor.tfs, median.cor.tfs.non, act.rep.thres.l, finalClassification = output.global.TFs, file = NULL)
    dev.off()

    
    ###############################
    # DESEq results, needed after #
    ###############################  
    # Run DESeq and valculate log2fc for genes, only needed here
    dd.filt <- DESeq(dd.filt)
    
    pdf(par.l$files_plotDiagnostic[3])

    #######################################
    # Correlation plots for the 3 classes #
    #######################################

    plot_classCorrelations(dd.filt, output.global.TFs, HOCOMOCO_mapping.df.overlap.exp, par.l, file = NULL)
    
    ###########################
    # DESEq2 diagnostic plots #
    ###########################
    plotDiagnosticPlots(dd.filt, conditionComparison, file = NULL, maxPairwiseComparisons = 0, alpha = 0.05, comparisonType = comparisonType, plotMA = TRUE)
    
    #############################
    # Density plots for each TF #
    #############################
    
    plot_density(t.cor.sel.matrix, t.cor.sel.matrix.non, file = NULL)
    
    dev.off()
    
} else {
    classesList.l = list(c())
}

output.global.TFs$yValue = transform_yValues(output.global.TFs$pvalueAdj, addPseudoCount = TRUE, nPermutations = par.l$nPermutations)
output.global.TFs.origReal = output.global.TFs

#########################################
# PLOT FOR DIFFERENT P VALUE THRESHOLDS #
#########################################

# Set the page dimensions to the maximum across all plotted variants
output.global.TFs.filteredSummary = dplyr::filter(output.global.TFs, pvalue <= max(par.l$significanceThresholds))

nTF_label = min(par.l$maxTFsToLabel, nrow(output.global.TFs.filteredSummary))


TFLabelSize = chooseTFLabelSize(nTF_label)

# Currently ignored, can be used to exchange weighted_meanDifference for something else
variableXAxis = "weighted_meanDifference"

thresholds = "NA" # not applicable here
if (par.l$plotRNASeqClassification) {
    thresholds  = par.l$allClassificationThresholds
}

allPlots.l = list()

flog.info(paste0("Looping through stringency thresholds and prepare plots..."))
for (thresholdCur in thresholds) {
    
    flog.info(paste0(" Stringency threshold for classification: ", thresholdCur))
    
    allPlots.l[[as.character(thresholdCur)]] = list()
    colnameClassificationFinal = paste0("classification_q", thresholdCur, "_final")
    
    for (significanceThresholdCur in par.l$significanceThresholds) {
        
        pValThrStr = as.character(significanceThresholdCur)
        
        for (showClasses in classesList.l) {
            
            output.global.TFs = output.global.TFs.origReal %>%
                mutate( pValueAdj_log10 = transform_yValues(pvalueAdj, addPseudoCount = TRUE, nPermutations = par.l$nPermutations),
                        pValue_log10    = transform_yValues(pvalue,    addPseudoCount = TRUE, nPermutations = par.l$nPermutations),
                        pValueAdj_sig   = pvalueAdj <= significanceThresholdCur,
                        pValue_sig      = pvalue <= significanceThresholdCur) 
            
            if (par.l$plotRNASeqClassification) {
                #output.global.TFs = filter(output.global.TFs, classification %in% showClasses)
                output.global.TFs = dplyr::filter(output.global.TFs, !!as.name(colnameClassificationFinal) %in% showClasses)
                
            }
            
            for (pValueStrCur in c("pvalue", "pvalueAdj")) {
                
                if (pValueStrCur == "pvalue") {
                    
                    pValueScoreCur = "pValue_log10"
                    pValueSigCur = "pValue_sig"
                    pValueStrLabel = "raw p-value"
                    
                    ggrepel_df = dplyr::filter(output.global.TFs, pValue_sig == TRUE)
                    maxPValue = max(output.global.TFs$pValue_log10, na.rm = TRUE)
                    
                } else {
                    
                    pValueScoreCur = "pValueAdj_log10"
                    pValueSigCur = "pValueAdj_sig"
                    pValueStrLabel = "adj. p-value"
                    
                    ggrepel_df = dplyr::filter(output.global.TFs, pValueAdj_sig == TRUE)
                    maxPValue = max(output.global.TFs$pValueAdj_log10, na.rm = TRUE)
                }
                
                # Increase the ymax a bit more
                ymax = max(transform_yValues(significanceThresholdCur, addPseudoCount = TRUE, nPermutations = par.l$nPermutations), maxPValue, na.rm = TRUE) * 1.1
                alphaValueNonSign = 0.3
                
                # Reverse here because negative values at left mean that the condition that has been specified in the beginning is higher. 
                # Reverse the rev() that was done before for this plot therefore to restore the original order
                labelsConditionsNew = rev(conditionComparison)
                
                g = ggplot()
                
                label_nameChange = 'TF activity higher in'
                if (comparisonMode != "pairwise") {
                    label_nameChange = 'Change with increasing\nvalues of conditionSummary'
                }
                
                if (par.l$plotRNASeqClassification) {

                    labelLegendCur = paste0("TF class (stringency: ", thresholdCur, ")")
                    g = g + geom_point(data = output.global.TFs, aes_string("weighted_meanDifference", pValueScoreCur, 
                                                                            alpha = pValueSigCur, size = "TFBS", fill = colnameClassificationFinal),  shape=21, stroke = 0.5, color = "black") +  
                        scale_fill_manual(labelLegendCur, values = par.l$colorCategories)
                    
                    g = g + geom_rect(aes(xmin = -Inf,xmax = 0,ymin = -Inf, ymax = Inf, color = par.l$colorConditions[2]),
                                      alpha = .3, fill = par.l$colorConditions[2], size = 0) +
                        geom_rect(aes(xmin = 0, xmax = Inf, ymin = -Inf,ymax = Inf, color = par.l$colorConditions[1]), alpha = .3, fill = par.l$colorConditions[1], size = 0) + 
                        scale_color_manual(name = label_nameChange, values = par.l$colorConditions, labels = conditionComparison)
                    
                } else {
                    
                    g = g + geom_point(data = output.global.TFs, aes_string("weighted_meanDifference", pValueScoreCur, alpha = pValueSigCur, size = "TFBS"), shape=21, stroke = 0.5, color = "black")
                    g = g + geom_rect(aes(xmin = -Inf, xmax = 0,   ymin = -Inf, ymax = Inf, fill = par.l$colorConditions[2]), alpha = .3) + 
                        geom_rect(aes(xmin = 0,    xmax = Inf, ymin = -Inf, ymax = Inf, fill = par.l$colorConditions[1]), alpha = .3)
                    g = g + scale_fill_manual(name = 'TF activity higher in', values = rev(par.l$colorConditions), labels = labelsConditionsNew)
                }
                
                g = g + ylim(-0.1,ymax) + 
                    ylab(paste0(transform_yValues_caption(), " (", pValueStrLabel, ")")) + 
                    xlab("weighted mean difference") + 
                    scale_alpha_manual(paste0(pValueStrLabel, " < ", significanceThresholdCur), values = c(alphaValueNonSign, 1), labels = c("no", "yes")) + 
                    geom_hline(yintercept = transform_yValues(significanceThresholdCur, addPseudoCount = TRUE, nPermutations = par.l$nPermutations), linetype = "dotted") 
                
                if (nrow(ggrepel_df) <= par.l$maxTFsToLabel) {
                    
                    if (par.l$plotRNASeqClassification) {
                        g = g +  geom_label_repel(data = ggrepel_df, aes_string("weighted_meanDifference", pValueScoreCur, label = "TF", fill = colnameClassificationFinal),
                                                  size = TFLabelSize, fontface = 'bold', color = 'white',
                                                  segment.size = 0.3, box.padding = unit(0.2, "lines"), max.iter = 5000,
                                                  label.padding = unit(0.2, "lines"), # how thick is connectin line
                                                  nudge_y = 0.05, nudge_x = 0,  # how far from center points
                                                  segment.alpha = .8, segment.color = par.l$plot_grayColor, show.legend = FALSE)
                    } else {
                        g = g +  geom_label_repel(data = ggrepel_df, aes_string("weighted_meanDifference", pValueScoreCur, label = "TF"),
                                                  size = TFLabelSize, fontface = 'bold', color = 'black',
                                                  segment.size = 0.3, box.padding = unit(0.2, "lines"), max.iter = 5000,
                                                  label.padding = unit(0.2, "lines"), # how thick is connectin line
                                                  nudge_y = 0.05, nudge_x = 0,  # how far from center points
                                                  segment.alpha = .8, segment.color = par.l$plot_grayColor, show.legend = FALSE)
                    }
                } else {
                    
                    flog.warn(paste0(" Not labeling significant TFs, maximum of ", par.l$maxTFsToLabel, " exceeded for ", pValThrStr, " and ", pValueStrCur))
                    
                    labelPlot = paste0("*TF labeling skipped because number of significant TFs\nexceeds the maximum of ", par.l$maxTFsToLabel, " (", nrow(ggrepel_df), ")")
                    
                    g = g + annotate("text", label = labelPlot, x = 0, y = ymax, size = 3)
                    
                } # end else
                
                g = g + theme_bw() + 
                    theme(axis.text.x = element_text(size=rel(1.5)), axis.text.y = element_text(size=rel(1.5)), 
                          axis.title.x = element_text(size=rel(1.5)), axis.title.y = element_text(size=rel(1.5)),
                          legend.title=element_text(size=rel(1.5)), legend.text=element_text(size=rel(1.5))) 
                
                if (par.l$plotRNASeqClassification) {
                    g = g + guides(alpha = guide_legend(override.aes = list(size=5), order = 2),
                                   fill = guide_legend(override.aes = list(size=5), order = 3),
                                   color = guide_legend(override.aes = list(size=5), order = 1))
                    
                    allPlots.l[[as.character(thresholdCur)]] [[pValThrStr]] [[paste0(showClasses,collapse = "-")]] [[pValueStrCur]] = g
                    
                } else {
                    g = g + guides(alpha = guide_legend(override.aes = list(size=5), order = 2),
                                   fill = guide_legend(override.aes = list(size=5), order = 3))
                    
                    allPlots.l[[as.character(thresholdCur)]] [[pValThrStr]] [[pValueStrCur]] = g
                }
                
            } # end separately for raw and adjusted p-values
        } # end for all showClasses
    } # end for different significance thresholds
    
    
    ####################
    # VOLCANO PLOT PDF #
    ####################
    height = width = max(nTF_label / 15 , par.l$volcanoPlot_minDimensions)
    
    if (par.l$plotRNASeqClassification) {
        
        pdf(file = par.l$file_plotVolcano, height = height, width = width, useDingbats = FALSE)
        plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
        message = paste0("This file is intentionally empty.\n",
                         "\nSee the identically named files with the suffix \".q\" for results.\n",
                         "\nIncreasing q-values indicate higher stringency\n and therefore more TFs are classified as undetermined.\n")
        text(x = 0.5, y = 0.5, message, cex = 0.9, col = "black")
        dev.off()
        
        fileCur = gsub(".pdf", paste0(".q", thresholdCur, ".pdf"), par.l$file_plotVolcano)
        pdf(file = fileCur, height = height, width = width, useDingbats = FALSE)
        
    } else { # Only one variant here
        pdf(file = par.l$file_plotVolcano, height = height, width = width, useDingbats = FALSE)
    }
    

    for (pValueStrCur in c("pvalueAdj", "pvalue")) {
        
        for (significanceThresholdCur in par.l$significanceThresholds) {
            
            for (showClasses in classesList.l) {
                
                if (par.l$plotRNASeqClassification) {
                    plot(allPlots.l[[as.character(thresholdCur)]] [[as.character(significanceThresholdCur)]] [[paste0(showClasses,collapse = "-")]] [[pValueStrCur]])
                } else {
                    plot(allPlots.l[[as.character(thresholdCur)]] [[as.character(significanceThresholdCur)]] [[pValueStrCur]])
                }
            } # end for each class
        } # end for each threshold
    } # end for both raw and adjusted p-values
    dev.off()
    
    
    
} # end for all stringency thresholds for the classification

#########################
# FINALLY, SAVE TO DISK #
#########################
output.global.TFs.origReal = dplyr::select(output.global.TFs.origReal, -one_of("permutation", "yValue"))

# TODO: Change 
# mutate_at(c("Diff_mean", "Diff_median", "Diff_mode", "Diff_skew", "Mean_l2FC", 
#"Median_l2FC", "Mode_l2FC", "skewness_l2FC"), signif, 3) %>%
 #   mutate_at(c("pvalue_raw", "adj_pvalue"), formatC, format = "g", digits = 3)

output.global.TFs.origReal.transf = dplyr::mutate_if(output.global.TFs.origReal, is.numeric, as.character)
write_tsv(output.global.TFs.origReal.transf, path = par.l$file_output_summary, col_names = TRUE)
saveRDS(allPlots.l, file = par.l$file_output_plots)

.printExecutionTime(start.time)
flog.info("Session info: ", sessionInfo(), capture = TRUE)
