start.time  <-  Sys.time()


#########################
# LIBRARY AND FUNCTIONS #
#########################


# Use the following line to load the Snakemake object to manually rerun this script (e.g., for debugging purposes)
# Replace {outputFolder} and {TF} correspondingly.
# snakemake = readRDS("{outputFolder}/LOGS_AND_BENCHMARKS/analyzeTF.{TF}.R.rds")

library("checkmate")
assertClass(snakemake, "Snakemake")
assertDirectoryExists(snakemake@config$par_general$dir_scripts)
source(paste0(snakemake@config$par_general$dir_scripts, "/functions.R"))

initFunctionsScript(packagesReq = NULL, minRVersion = "3.1.0", warningsLevel = 1, disableScientificNotation = TRUE)
checkAndLoadPackages(c("tidyverse", "futile.logger", "DESeq2", "vsn", "modeest", "checkmate", "limma", "geneplotter", "RColorBrewer", "tools"), verbose = FALSE)

########################################################################
# SAVE SNAKEMAKE S4 OBJECT THAT IS PASSED ALONG FOR DEBUGGING PURPOSES #
########################################################################

# Use the following line to load the Snakemake object to manually rerun this script (e.g., for debugging purposes)
# Replace {outputFolder} and {TF} correspondingly.
# snakemake = readRDS("{outputFolder}/LOGS_AND_BENCHMARKS/3.analyzeTF.{TF}.R.rds")

createDebugFile(snakemake)

###################
#### PARAMETERS ###
###################

par.l = list()

par.l$verbose = TRUE
par.l$log_minlevel = "INFO"
par.l$maxPairwiseComparisonsDiagnosticPermutations = 2
par.l$pseudocountAddition = 1
par.l$minNoDatapoints = 10
par.l$roundLog2FCDigits = 5

#####################
# VERIFY PARAMETERS #
#####################

assertClass(snakemake, "Snakemake")

## INPUT ##
assertList(snakemake@input, min.len = 1)
assertSubset(names(snakemake@input), c("", "overlapFile", "sampleDataR", "peakFile", "peakFile2", "normFacs", "condComp"))

par.l$file_input_peakTFOverlaps  = snakemake@input$overlapFile
assertFileExists(par.l$file_input_peakTFOverlaps, access = "r")

par.l$file_input_metadata = snakemake@input$sampleDataR
assertFileExists(par.l$file_input_metadata, access = "r")

par.l$file_input_peaks = snakemake@input$peakFile
assertFileExists(par.l$file_input_peaks, access = "r")

par.l$file_input_peak2 = snakemake@input$peakFile2
assertFileExists(par.l$file_input_peak2, access = "r")

par.l$file_input_normFacs = snakemake@input$normFacs
assertFileExists(par.l$file_input_normFacs, access = "r")

par.l$file_input_conditionComparison = snakemake@input$condComp
assertFileExists(par.l$file_input_conditionComparison, access = "r")


## OUTPUT ##
assertList(snakemake@output, min.len = 1)
assertSubset(names(snakemake@output), c("", "outputTSV", "outputPermTSV", "outputRDS", "plot_diagnostic"))

par.l$file_output_summaryAll      = snakemake@output$outputTSV
par.l$file_outputPerm_summaryAll  = snakemake@output$outputPermTSV
par.l$file_output_summaryStats    = snakemake@output$outputRDS
par.l$file_output_plot_diagnostic = snakemake@output$plot_diagnostic



## WILDCARDS ##
assertList(snakemake@wildcards, min.len = 1)
assertSubset(names(snakemake@wildcards), c("", "TF"))

par.l$TF = snakemake@wildcards$TF
assertCharacter(par.l$TF, len = 1, min.chars = 1)

## CONFIG ##
assertList(snakemake@config, min.len = 1)

par.l$designFormula = snakemake@config$par_general$designContrast
checkAndLogWarningsAndErrors(par.l$designFormula, checkCharacter(par.l$designFormula, len = 1, min.chars = 3))

par.l$designFormulaVariableTypes = snakemake@config$par_general$designVariableTypes
checkAndLogWarningsAndErrors(par.l$designFormulaVariableTypes, checkCharacter(par.l$designFormulaVariableTypes, len = 1, min.chars = 3))

par.l$nPermutations = snakemake@config$par_general$nPermutations
assertIntegerish(par.l$nPermutations, lower = 0)

par.l$conditionComparison  = snakemake@config$par_general$conditionComparison
checkAndLogWarningsAndErrors(par.l$conditionComparison, checkCharacter(par.l$conditionComparison, len = 1))


par.l$debugMode = setDebugMode(snakemake@config$par_general$debugMode)


## PARAMS ##
assertList(snakemake@params, min.len = 1)
assertSubset(names(snakemake@params), c("", "doCyclicLoess", "allBAMS", "debugFile"))

par.l$doCyclicLoess = as.logical(snakemake@params$doCyclicLoess)
assertFlag(par.l$doCyclicLoess)

par.l$allBAMS = snakemake@params$allBAMS
for (fileCur in par.l$allBAMS) {
  checkAndLogWarningsAndErrors(fileCur, checkFileExists(fileCur))
}

## LOG ##
assertList(snakemake@log, min.len = 1)
par.l$file_log = snakemake@log[[1]]


allDirs = c(dirname(par.l$file_output_summaryAll), 
            dirname(par.l$file_outputPerm_summaryAll),
            dirname(par.l$file_output_summaryStats), 
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
    sampleData.l = readRDS(par.l$file_input_metadata)
    normFacs = readRDS(par.l$file_input_normFacs)
    
    nTFBS = length(readLines(par.l$file_input_peakTFOverlaps))
    if (nTFBS > 0) {
        overlapsAll.df = read_tidyverse_wrapper(par.l$file_input_peakTFOverlaps, type = "tsv", col_names = TRUE, col_types = cols(), comment = "#")
    }
    peaks.df = read_tidyverse_wrapper(par.l$file_input_peak2, type = "tsv", col_types = cols())
    peaksFiltered.df = readRDS(par.l$file_input_peaks)
    conditionComparison = readRDS(par.l$file_input_conditionComparison)
    
    save(list = ls(), file = snakemake@params$debugFile)
    flog.info(paste0("File ", snakemake@params$debugFile, " has been saved. You may use it for trouble-shooting and debugging, see the Documentation for more details."))
    
}

#################
# READ METADATA #
#################

sampleData.l = readRDS(par.l$file_input_metadata)

if (length(sampleData.l) == 0) {
    message = "Length of sampleData.l list is 0 but is has to be at least 1. Rerun the rule DiffPeaks."
    checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
}


# Adjust the number of permutations in case less have been computed
if (par.l$nPermutations + 1 < length(sampleData.l)) {
  message = paste0("In the output objects, more permutations seem to be stored. They will be ignored and the currently specified value of nPermutations will be used")
  checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
} else if (par.l$nPermutations + 1 > length(sampleData.l)) {
  valueNew = length(sampleData.l) - 1
  message = paste0("The value of the parameter nPermutations differs from what is saved in the output objects. The value of nPermutations will be adjusted to ", valueNew)
  checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
  par.l$nPermutations = valueNew
}


sampleData.df = sampleData.l[["permutation0"]]
colnamesNew = sampleData.df$SampleID[which(sampleData.df$bamReads %in% par.l$allBAMS)]
if (length(colnamesNew) != nrow(sampleData.df)) {
    message = "Could not grep sampleIDs from filenames."
    checkAndLogWarningsAndErrors(NULL,  message, isWarning = FALSE)
}


designComponents.l = checkDesignIntegrity(snakemake, par.l, sampleData.df)

components3types   = designComponents.l$types
variableToPermute  = designComponents.l$variableToPermute

# Which of the two modes should be done, pairwise or quantitative?
comparisonMode = "quantitative"
if (components3types["conditionSummary"] == "logical" | components3types["conditionSummary"] == "factor") {
    comparisonMode = "pairwise"
}


# Initiate data structures that are populated hereafter
TF_output.df  = tribble(~permutation, ~TF, ~chr, ~MSS, ~MES, ~TFBSID, ~strand, ~peakID, ~limma_avgExpr, ~l2FC, ~limma_B, ~limma_t_stat, ~DESeq_ldcSE, ~DESeq_stat, ~DESeq_baseMean, ~pval, ~pval_adj)


outputSummary.df  = tribble(~permutation, ~TF, ~Pos_l2FC, ~Mean_l2FC, ~Median_l2FC, ~Mode_l2FC, ~sd_l2FC, ~pvalue_raw, ~skewness_l2FC, ~T_statistic, ~TFBS_num)


#####################
# READ OVERLAP FILE #
#####################

# Check number of lines. If file is empty, the TF has to be skipped
nTFBS = length(readLines(par.l$file_input_peakTFOverlaps))
if (nTFBS > 0) {
    
    overlapsAll.df = read_tidyverse_wrapper(par.l$file_input_peakTFOverlaps, type = "tsv",
                                            col_names = TRUE, col_types = cols(), comment = "#")
    nTFBS = nrow(overlapsAll.df)
    
    colnames(overlapsAll.df) = c("annotation", "chr","MSS","MES", "strand","length", colnamesNew)
    
    overlapsAll.df = overlapsAll.df %>%
        dplyr::mutate(TFBSID = paste0(chr,":", MSS, "-",MES),
                      mean = apply(dplyr::select(overlapsAll.df, one_of(colnamesNew)), 1, mean), 
                      peakID = sapply(strsplit(overlapsAll.df$annotation, split = "_", fixed = TRUE),"[[", 1)) %>%
        dplyr::distinct(TFBSID, .keep_all = TRUE) %>%
        dplyr::select(-one_of("length"))
    
    skipTF = FALSE
} else {
  skipTF = TRUE
}


if (nTFBS >= par.l$minNoDatapoints) {

    # Create formula based on user-defined design
    designFormula = convertToFormula(par.l$designFormula, colnames(sampleData.df))
    formulaVariables = attr(terms(designFormula), "term.labels")
    # Extract the variable that defines the contrast. Always the last element in the formula
    variableToPermute = formulaVariables[length(formulaVariables)]
    
  # Group by peak ID: To avoid biases and dependencies based on TFBS clustering within peaks, we then select the TFBS per TF per peak with the highest average read count across all samples.
  coverageAll_grouped.df = overlapsAll.df %>%
    dplyr::group_by(peakID) %>%
    dplyr::slice(which.max(mean)) %>%
    dplyr::ungroup()
  
  TF.table.m = as.matrix(coverageAll_grouped.df[,sampleData.df$SampleID])
  colnames(TF.table.m) = sampleData.df$SampleID
  rownames(TF.table.m) = coverageAll_grouped.df$TFBSID
  
} else {
  message <- paste0("The number of overlaps (TFBS) with the peaks is smaller than ", par.l$minNoDatapoints, ". This TF will be ignored in subsequent steps")
  checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
  skipTF = TRUE
}

if (skipTF) {
  
  # write a dummy pdf file
  pdf(par.l$file_output_plot_diagnostic)
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  message = paste0("Insufficient data to run analysis.\n", message, "\nThis TF will be ignored in subsequent steps.")
  text(x = 0.5, y = 0.5, message, cex = 1.6, col = "red")
  dev.off()
  
  
  TF_outputInclPerm.df = as.data.frame(matrix(nrow = 0, ncol = 2 + par.l$nPermutations + 1))
  colnames(TF_outputInclPerm.df) = c("TF", "TFBSID", paste0("log2fc_perm", 0:par.l$nPermutations))
  
  
} else {

  normFacs = readRDS(par.l$file_input_normFacs)
  
  
  TF.cds = tryCatch( {
      
      # create Deseq object from the TF specific data
      # The correct order is already enforced due to the creation of the TF.table.m matrix before that is sorted after sampleData.df
      TF.cds <- DESeqDataSetFromMatrix(countData = TF.table.m,
                                       colData = sampleData.df,
                                       design = designFormula)
      
      # Recent versions of DeSeq seem to do this automatically, whereas older versions don't, so enforce it here
      if (!identical(colnames(TF.cds), colnames(TF.table.m))) {
          colnames(TF.cds) = colnames(TF.table.m)
      }
      if (!identical(rownames(TF.cds), rownames(TF.table.m))) {
          rownames(TF.cds) = rownames(TF.table.m)
      }
      
      TF.cds
      
  }, error = function(e) {
      errorMessage = "Could not initiate DESeq."
      checkAndLogWarningsAndErrors(NULL, errorMessage, isWarning = FALSE)
  }
  )
  
  
  # We here again take the gene-specific normalization factors from the previously calculated ones based on the peaks,
  # but subset this to only those corresponding to the TF of interest
  
  
  # Sanity check
  if (length(which(!coverageAll_grouped.df$peakID %in% rownames(normFacs))) > 0) {
      errorMessage = "Inconsistency detected between the normalization factor rownames and the object coverageAll_grouped.df"
      checkAndLogWarningsAndErrors(NULL,  errorMessage, isWarning = FALSE)
  }
  
  if (!identical(colnames(TF.cds), colnames(normFacs))) {
      errorMessage = "Column names different between TF.cds and normFacs"
      checkAndLogWarningsAndErrors(NULL,  errorMessage, isWarning = FALSE)
  }
  
  if (!identical(coverageAll_grouped.df$TFBSID, rownames(TF.cds))) {
      errorMessage = "Row names different between coverageAll_grouped.df and TF.cds"
      checkAndLogWarningsAndErrors(NULL,  errorMessage, isWarning = FALSE)
  }
  
  
  if (par.l$doCyclicLoess) {
      
      # coverageAll_grouped.df$TFBSID and rownames(TF.cds) are identical
      matchingOrder = match(coverageAll_grouped.df$peakID,rownames(normFacs))
      normalizationFactors(TF.cds) <- normFacs[matchingOrder,, drop = FALSE]
      
  } else {
      sizeFactors(TF.cds) = normFacs
  }
  
  # low RC, check by rowMean
  TF.cds.filt = TF.cds[rowMeans(DESeq2::counts(TF.cds)) > 0, ]
  
  nPeaks = nrow(TF.cds.filt)
  
  # Preallocate data frame so no expensive reallocation has to be done
  log2fc.m = matrix(NA, nrow = nPeaks , ncol = par.l$nPermutations + 1)

  peaks.df = read_tidyverse_wrapper(par.l$file_input_peak2, type = "tsv", col_types = cols())
  
  peaksFiltered.df = readRDS(par.l$file_input_peaks)
  
  # READ FILE HERE
  conditionComparison = readRDS(par.l$file_input_conditionComparison)
  
  ################################
  # ITERATE THROUGH PERMUTATIONS #
  ################################

    # Calculate the log2 counts once
    if (par.l$nPermutations > 0) {
      
      # Generate normalized counts for limma analysis
      countsRaw         = DESeq2::counts(TF.cds.filt, norm = FALSE)
      countsNorm        = DESeq2::counts(TF.cds.filt, norm = TRUE)
      countsNorm.transf = log2(countsNorm + par.l$pseudocountAddition)
      rownames(countsNorm.transf) = rownames(TF.cds.filt)
      
    }
    
    # dplyr::mutate(TF = par.l$TF) gives the following weird error message: Error: Unsupported type NILSXP for column "TF"
    TFCur = par.l$TF
    
    for (permutationCur in 0:par.l$nPermutations) {
      
        if (permutationCur > 0 & (permutationCur %% 10 == 0 | permutationCur == par.l$nPermutations)) {
            flog.info(paste0("Running permutation ", permutationCur))
        } else {
            flog.info(paste0("Running for real data "))
        }
        
      sampleData.df = sampleData.l[[paste0("permutation", permutationCur)]]
      
      ##############################
      # RUN EITHER LIMMA OR DESEQ2 #
      ##############################
      if (par.l$nPermutations > 0) {
        
          # model.matrix uses the first level in the specified column as reference, and so the corresponding column name and values are relative to that reference.
          # That is, if the levels are "GMP" and "MPP", then all log2 fc will be the log2fc of MPP as compared to GMP.
          # whatever comes first for model.matrix is taken as first value, then log2fc is of the second condition over the first
        fit <- eBayes(lmFit(countsNorm.transf, design = model.matrix(designFormula, data = sampleData.df)))
        results.df <- topTable(fit, coef = colnames(fit$design)[ncol(fit$design)], number = Inf, sort.by = "none")
        
        final.TF.df = tibble("TFBSID"      = rownames(results.df), 
                                 "limma_avgExpr"     = results.df$AveExpr,
                                 "l2FC"        = results.df$logFC,
                                 "limma_B"           = results.df$B,
                                 "limma_t_stat"      = results.df$t,
                                 "pval"        = results.df$P.Value, 
                                 "pval_adj"    = results.df$adj.P.Val,
                                 "DESeq_ldcSE" = NA, 
                                 "DESeq_stat" = NA,
                                 "DESeq_baseMean" = NA
        )
        
      } else {
        
        sampleData.df = sampleData.l[[paste0("permutation0")]]
        
        # We already set the factors for conditionSummary explicitly. The reference level is the first level for DeSeq. 
        # Run the local fit first, if that throws an error try the default fit type
        
        TF.cds.filt = tryCatch( {
          suppressMessages(DESeq(TF.cds.filt,fitType = 'local'))
          
        }, error = function(e) {
          message = "Could not run DESeq with local fitting, retry with default fitting type..."
          checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
          
          TF.cds.filt = tryCatch( {
            suppressMessages(DESeq(TF.cds.filt))
            
          }, error = function(e) {
            errorMessage <<- "Could not run DESeq with regular fitting either, set all values to NA."
            checkAndLogWarningsAndErrors(NULL, errorMessage, isWarning = TRUE)
          }
          )
          
          TF.cds.filt
          
        }
        )
        
        if (class(TF.cds.filt) == "character") {
          
          skipTF = TRUE
          TF_outputInclPerm.df = as.data.frame(matrix(nrow = 0, ncol = 2 + par.l$nPermutations + 1))
          colnames(TF_outputInclPerm.df) = c("TF", "TFBSID", paste0("log2fc_perm", 0:par.l$nPermutations))
        }
        
        if (!skipTF) {
            
            #Enforce the correct order of the comparison
            if (comparisonMode == "pairwise") { 
                
                contrast = c(variableToPermute, conditionComparison[1], conditionComparison[2])
                
            } else {
                
                # Same as without specifying contrast at all
                contrast = list(variableToPermute)
                

            }
            
            res_DESeq = DESeq2::results(TF.cds.filt, contrast = contrast)
            res_DESeq.df <- as.data.frame(res_DESeq)
           
            
            final.TF.df = tibble("TFBSID"    = rownames(res_DESeq.df), 
                                     "DESeq_baseMean" = res_DESeq.df$baseMean,
                                     "l2FC"     = res_DESeq.df$log2FoldChange,
                                     "DESeq_ldcSE"    = res_DESeq.df$lfcSE,
                                     "DESeq_stat"     = res_DESeq.df$stat,
                                     "pval"     = res_DESeq.df$pvalue, 
                                     "pval_adj" = res_DESeq.df$padj,
                                     "limma_avgExpr" = NA,
                                     "limma_B" = NA,
                                     "limma_t_stat" = NA
            )
        }
      }
      
      ##################################
      # SUMMARIZE AND MERGE WITH PEAKS #
      ##################################
      if (!skipTF) {
          order  = c("permutation", "TF", "chr", "MSS", "MES", "TFBSID", "strand", "peakID", "l2FC", "limma_avgExpr", "limma_B", "limma_t_stat", "DESeq_ldcSE", "DESeq_stat", "DESeq_baseMean", "pval", "pval_adj")
     
          TF_outputCur.df = final.TF.df %>%
            dplyr::left_join(overlapsAll.df,by = c("TFBSID")) %>%
            dplyr::left_join(peaksFiltered.df, by = c("peakID" = "annotation")) %>%
            dplyr::rename(chr = chr.x) %>%
            dplyr::filter(!(is.na(limma_avgExpr) & is.na(DESeq_baseMean))) %>%
            dplyr::arrange(chr)  %>%
            dplyr::mutate(TF = TFCur, permutation = permutationCur, l2FC = signif(l2FC, par.l$roundLog2FCDigits)) %>%
            dplyr::select(one_of(order)) 
          
    
          if (permutationCur == 0) {
            
            TF_output.df = rbind(TF_output.df, TF_outputCur.df)
            
            ######################
            ## DIAGNOSTIC PLOTS ##
            ######################
  
            pdf(par.l$file_output_plot_diagnostic)
            if (par.l$nPermutations == 0) {
              plotDiagnosticPlots(TF.cds.filt, conditionComparison, file = NULL, maxPairwiseComparisons = 0,  plotMA = FALSE) 
            } else {
                
              plotDiagnosticPlots(fit, conditionComparison, file = NULL, maxPairwiseComparisons = 0, plotMA = FALSE, counts.raw = countsRaw, counts.norm = countsNorm) 
            }
            
            
            
            # Density plot
            conditionComparison = paste0(conditionComparison[1], " vs. ", conditionComparison[2])
            xlabLabel = paste0(" log2 FC ", conditionComparison)
            TF_dens = ggplot(peaks.df, aes(l2FC)) + geom_density(aes(fill = "Peaks" ),alpha = .5, color = "black") +
              geom_density(data = final.TF.df, aes(x = l2FC, fill = "TF"),size = 1, alpha = .7) +
              xlab(xlabLabel) +
              theme(axis.text.x = element_text(face = "bold", color = "black", size = 12),
                    axis.text.y = element_text(face = "bold", color = "black", size = 12),
                    axis.title.x = element_text(face = "bold", colour = "black", size = 10, margin = margin(25,0,0,0)),
                    axis.title.y = element_text(face = "bold", colour = "black", size = 10, margin = margin(0,25,0,0)),
                    axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank(),
                    legend.position = c(0.9,0.9),
                    legend.justification = "center",
                    legend.title = element_blank()) +
              scale_fill_manual(values = c("Peaks" = "grey50" , "TF" = "blue"), labels = c("PEAKS", par.l$TF))
            
            plot(TF_dens)
            
            # create ecdf plots for each TF
            ECDF_TF = ggplot(final.TF.df, aes(l2FC)) +
              stat_ecdf(aes(colour = "TF")) +
              stat_ecdf(data = peaks.df, aes(x = l2FC, colour = "Peaks")) +
              xlab(xlabLabel) +
              scale_fill_manual(values = c("Peaks" = "grey50" , "TF" = "blue"), labels = c("PEAKS", par.l$TF))
            
            plot(ECDF_TF)
            
            dev.off()
            
          } 
          
          # Execute this ALWAYS, as also the values for the real data should be stored for easier later retrieval
    
          log2fc.m[,permutationCur + 1] = TF_outputCur.df$l2FC
    
          # d) Comparisons between peaks and binding sites
          
          # Check the version of modeest, because version 2.3.2 introduced an implementation change that breaks things
          
          if (nrow(dplyr::filter(final.TF.df, !is.na(l2FC))) > 0) {
              
              if (packageVersion("modeest") < "2.3.2") {
                  
                  modeNum     = mlv(final.TF.df$l2FC, method = "mfv", na.rm = TRUE)
                  
                  stopifnot(is.list(modeNum))
                  l2fc_mode = ifelse(is.null(modeNum$M), NA, modeNum$M)
                  l2fc_skewness = ifelse(is.null(modeNum$skewness), NA, modeNum$skewness)
              } else {
                  
                  l2fc_mode = mlv(final.TF.df$l2FC, method = "mfv", na.rm = TRUE)[1]
                  l2fc_skewness = skewness(final.TF.df$l2FC, na.rm = TRUE)[1]
              }
              
          } else {
              
              message = paste0("l2fc values could not calculated for any TFBS, this may happen for complex design formulas in combination with particular permutations. For permutation ", permutationCur, ", all values are NA.")
              checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
              # For rare occassions, if only NA data are available, set l2fc_mode and l2fc_skewness to NA also
              l2fc_mode = l2fc_skewness = NA
          }
          
         
          
          
          # We have to filter now by NAs because for some permutations, limma might have been unable to calculate the coefficients
          final.TF.filtered.df = dplyr::filter(final.TF.df, !is.na(l2FC))
          nRowFilt = nrow(final.TF.filtered.df)
          if (nRowFilt > 0) {
            
            if (nRowFilt >= par.l$minNoDatapoints) {
              Ttest   = t.test(final.TF.df$l2FC, peaks.df$l2FC)
              tTest_pVal = Ttest$p.value
              tTest_stat = Ttest$statistic[[1]]
            } else {
              message = paste0("Skipping t-test due to insufficient number of TFBS (<", par.l$minNoDatapoints, "). The columns T_statistic and pvalue_raw will be set to NA.")
              checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
              tTest_pVal = tTest_stat = NA
            }
            
            
            outputSummary.df = add_row(outputSummary.df,
                                       permutation     = permutationCur,
                                       TF              = par.l$TF,
                                       Pos_l2FC        = nrow(final.TF.df[final.TF.df$l2FC > 0,]) / nrow(final.TF.df),
                                       Mean_l2FC       = mean(final.TF.df$l2FC, na.rm = TRUE),
                                       Median_l2FC     = median(final.TF.df$l2FC, na.rm = TRUE),
                                       Mode_l2FC       = l2fc_mode,
                                       sd_l2FC         = sd(final.TF.df$l2FC, na.rm = TRUE),
                                       pvalue_raw      = tTest_pVal,
                                       skewness_l2FC   = l2fc_skewness, 
                                       T_statistic     = tTest_stat, 
                                       TFBS_num        = nrow(final.TF.df)
            )
            
          } else {
            
            # Rest will be NA
            outputSummary.df = add_row(outputSummary.df,
                                       permutation     = permutationCur,
                                       TF              = par.l$TF,
                                       TFBS_num        = nrow(final.TF.df)               
            )
          }
          
      } else { # if skipTF
          
          message <- paste0("DESeq could not be run.")
          pdf(par.l$file_output_plot_diagnostic)
          plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
          message = paste0("Insufficient data to run analysis.\n", message, "\nThis TF will be ignored in subsequent steps.")
          text(x = 0.5, y = 0.5, message, cex = 1.6, col = "red")
          dev.off()
          
      }
    }  # end for each permutation
  
  if (!skipTF) {
      TF_outputInclPerm.df = as_tibble(log2fc.m) %>%
          add_column(TF = TFCur , TFBSID = TF_outputCur.df$TFBSID, .before = 1)
      
      colnames(TF_outputInclPerm.df)[3:ncol(TF_outputInclPerm.df)] = paste0("log2fc_perm", 0:par.l$nPermutations)
  } 
 

} # end if !skipTF

# Do it separately for each column because different rounding schemes might be needed
TF_output.df = TF_output.df %>% 
    mutate_at(c("l2FC", "limma_avgExpr", "limma_B", "limma_t_stat", "DESeq_ldcSE", "DESeq_stat", "DESeq_baseMean"), signif, 3) %>%
    mutate_at(c("pval", "pval_adj"), formatC, format = "g", digits = 3)


# TF_output.df = mutate_if(TF_output.df, is.numeric, as.character)
write_tsv(TF_output.df,     path = par.l$file_output_summaryAll)

# All numeric columns can be treated the same way
TF_outputInclPerm.df = mutate_if(TF_outputInclPerm.df, is.numeric, signif, 2)
write_tsv(TF_outputInclPerm.df, path = par.l$file_outputPerm_summaryAll, col_names = FALSE)
saveRDS(outputSummary.df,   file = par.l$file_output_summaryStats)

# saveRDS(res_DESeq.l, file = par.l$file_output_DESeq)

.printExecutionTime(start.time)

flog.info("Session info: ", sessionInfo(), capture = TRUE)
