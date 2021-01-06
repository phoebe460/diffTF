
initFunctionsScript <- function(packagesReq = NULL, minRVersion = "3.1.0", warningsLevel = 1, disableScientificNotation = TRUE, verbose = FALSE) {
    
    checkAndLoadPackages(c("checkmate"), verbose = verbose)
    assert(checkNull(packagesReq), checkCharacter(packagesReq, min.len = 1, min.chars = 1))
    assertCharacter(minRVersion, len = 1)
    assertInt(warningsLevel, lower = 0, upper = 2)
    assertFlag(disableScientificNotation)
    assertFlag(verbose)
    
    clearOpenDevices()
    
    
    # No annoying strings as factors by default
    options(stringsAsFactors = FALSE)
    
    # Print warnings as they occur
    options(warn = warningsLevel)
    
    # Just print 200 lines instead of 99999
    options(max.print = 200)
    
    # Disable scientific notation
    if (disableScientificNotation) options(scipen = 999)
    
    
    # We need at least R version 3.1.0 to continue
    stopifnot(getRversion() >= minRVersion)
    
    if (is.null(packagesReq)) {
        #packagesReq = .loadAllLibraries()
    }
    
    
    .detachAllPackages()
    
    checkAndLoadPackages(packagesReq, verbose = verbose)
    
}

checkAndLoadPackages <- function(packages, verbose = FALSE) {
    
    .checkAndInstallMissingPackages(packages, verbose = verbose)
    
    for (packageCur in packages) {
        suppressMessages(library(packageCur, character.only = TRUE))
    }
    
}


startLogger <- function(logfile, level, removeOldLog = TRUE, appenderName = "consoleAndFile", verbose = TRUE) {
    
    checkAndLoadPackages(c("futile.logger"), verbose = FALSE)
    
    assertSubset(level, c("TRACE", "DEBUG", "INFO", "WARN", "ERROR", "FATAL"))
    assertFlag(removeOldLog)
    assertSubset(appenderName, c("console", "file", "consoleAndFile"))
    assertFlag(verbose)
    
    if (appenderName != "console") {
        assertDirectory(dirname(logfile), access = "w")
        if (file.exists(logfile)) {
            file.remove(logfile)
        }
    }
    
    # LEVELS: TRACE, DEBUG, INFO, WARN, ERROR, FATAL
    invisible(flog.threshold(level))
    
    
    if (appenderName == "console") {
        invisible(flog.appender(appender.console()))
    } else if (appenderName == "file")  {
        invisible(flog.appender(appender.file(file = logfile)))
    } else {
        invisible(flog.appender(appender.tee(file = logfile)))
    }
    
    
}

stopifnot_custom <- function(testForTruth, directory) {
    
    assert_flag(testForTruth)
    if (!testForTruth) {
        
        flog.error()
    }
}

printParametersLog <- function(par.l, verbose = FALSE) {
    
    checkAndLoadPackages(c("futile.logger"), verbose = verbose)  
    assertList(par.l)
    flog.info(paste0("PARAMETERS:"))
    for (parCur in names(par.l)) {
        
        flog.info(paste0(" ", parCur, "=",  paste0(par.l[[parCur]], collapse = ",")))
        
    }
}

read_tidyverse_wrapper <- function(file, type = "tsv", ncolExpected = NULL, minRows = 0, ...) {
  
  assertSubset(type, c("csv", "csv2", "tsv", "delim"))
  
  start = Sys.time()
  flog.info(paste0("Reading file ", file))
  
  
  if (type == "tsv") {
    tbl = read_tsv(file, ...)
  } else if (type == "csv") {
    tbl = read_csv(file, ...)
  } else if (type == "csv2") {
    tbl = read_csv2(file, ...)
  } else if (type == "delim") {
    tbl = read_delim(file, ...)
  }
  
 
  if (nrow(problems(tbl)) > 0) {
    flog.fatal(paste0("Parsing errors: "), problems(tbl), capture = TRUE)
    stop("Error when parsing the file ", file, ", see errors above")
  }
  
  
  if (nrow(tbl) == 0) {
    message = paste0("The file ", file, " is unexpectedly empty.")
    checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
  }
  
  if (!is.null(ncolExpected)) {
    if (! ncol(tbl) %in% ncolExpected) {
      message = paste0("The file ", file, " does not have the expected number of ", ncolExpected, " columns, but instead ", ncol(tbl), ".")
      checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
    }
  }
  
  if (minRows > 0) {
      if (nrow(tbl) < minRows) {
          message = paste0("The file ", file, " does not have the expected minimum number of rows. Expected at least ", minRows, ", but found only ",nrow(tbl), ".")
          checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
      }
  }
  
  
  .printExecutionTime(start)
  tbl
}


###########################################
# PACKAGE LOADING AND DETACHING FUNCTIONS #
###########################################

.checkAndInstallMissingPackages <- function(packages.vec, verbose = FALSE) {
    
    if (verbose) cat("Trying to automatically install missing packages. If this fails, install them manually...\n")
    
    packagesToInstall = setdiff(packages.vec, rownames(installed.packages()))
    
    
    if (length(packagesToInstall) > 0) {
        cat("Could not find the following packages: ", paste( packagesToInstall , collapse = ", "), "\n")
        install.packages(packagesToInstall, repos = "http://cran.rstudio.com/")  
        
        if (getRversion() < "3.5.0") { 
            
            source("http://bioconductor.org/biocLite.R")
            for (packageCur in packagesToInstall) {
                biocLite(packageCur, suppressUpdates = TRUE)
            }
        } else {
            
            if (!requireNamespace("BiocManager", quietly = TRUE))
                install.packages("BiocManager")
            
            for (packageCur in packagesToInstall) {
                BiocManager::install(packageCur, update = FALSE, ask = FALSE)
            }
        }
        
        
        
        
        
    } else {
        if (verbose) cat("All packages are already installed\n")
    }
    
}


.detachAllPackages <- function() {
    
    basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
    
    package.list <- search()[ifelse(unlist(gregexpr("package:",search())) == 1,TRUE,FALSE)]
    
    package.list <- setdiff(package.list,basic.packages)
    
    if (length(package.list) > 0)  for (package in package.list) detach(package, character.only = TRUE)
    
}



clearOpenDevices <- function() {
    
    while (length(dev.list()) > 0) {
        dev.off()
    }
}


createFileList <- function(directory, pattern, recursive = FALSE, ignoreCase = FALSE, verbose = TRUE) {
    
    assertCharacter(directory, min.chars = 1, any.missing = FALSE, len = 1)
    assertCharacter(pattern, min.chars = 1, any.missing = FALSE, len = 1)
    assertFlag(recursive)  
    assertFlag(ignoreCase)  
    assertFlag(verbose)    
    
    assertDirectoryExists(directory)
    
    # Multiple patterns are now supported, integrate over them
    patternAll = strsplit(pattern, ",")[[1]]
    assertCharacter(patternAll, min.len = 1)
    
    if (verbose) cat("Found ", length(patternAll), " distinct pattern(s) in pattern string.\n")
    
    nFilesToProcessTotal = 0
    filesToProcess.vec = c()
    
    for (patternCur in patternAll) {
        
        # Replace wildcards by functioning patterns (such as .)
        patternMod = glob2rx(patternCur)
        
        # Remove anchoring at beginning and end
        patternMod = substr(patternMod, 2, nchar(patternMod) - 1)
        
        filesToProcessCur.vec = list.files(path = directory, pattern = patternMod, full.names = TRUE, recursive = recursive, ignore.case = ignoreCase)
        filesToProcess.vec = c(filesToProcess.vec, filesToProcessCur.vec)
        
        if (verbose) cat("Search for files with pattern \"", patternCur, "\" in directory ", directory, " (case insensitive:", ignoreCase, ")\n", sep = "")
        
        nFilesToProcessTotal = nFilesToProcessTotal + length(filesToProcessCur.vec)
    }
    
    
    
    if (nFilesToProcessTotal == 0) {
        stop(paste0("No files to process in folder ", directory, " that fulfill the desired criteria (", patternCur, ")."))
    } else {
        
        if (verbose) cat("The following", nFilesToProcessTotal, "files were found:\n", paste0(filesToProcess.vec, collapse = "\n "))
    }
    
    if (verbose) cat("\n")
    
    filesToProcess.vec
    
}

getComparisonFromDeSeqObject <- function(DeSeq.obj, designFormula, datatypeVariableToPermute = "factor") {
    
    assertClass(DeSeq.obj, "DESeqDataSet")
    
    mcol.df = mcols(results(DeSeq.obj))
    resultsRow = mcol.df[which(mcol.df$type == "results"),][1,"description"]
    split1 = strsplit(resultsRow, split = ":", fixed = TRUE)[[1]][2]
    
    if (datatypeVariableToPermute == "factor") {
        
        formulaElements = gsub(pattern = "[~+*]", replacement = "", x = designFormula)
        elems = strsplit(trimws(formulaElements), split = "\\s+", perl = TRUE)[[1]]
        
        # Remove names of variables from split1 variable
        split1New = split1
        for (varCur in elems) {
            split1New = trimws(gsub(pattern = varCur, replacement = "", split1New))
        }
        
        splitNew = trimws(gsub(pattern = "vs", replacement = "", split1New))
        splitFinal = strsplit(splitNew, split = "\\s+", perl = TRUE)[[1]]
        
        
    } else if (datatypeVariableToPermute == "integer") {
        
        # Here, we manually assign neg. and pos. change
        splitFinal = trimws(split1)
        splitFinal = paste0(splitFinal, " (", c("pos.", "neg."), " change)")
        
    }
    
    assertVector(splitFinal, len = 2)
    return(splitFinal)
}


permuteSampleTable <- function(file_sampleTable, conditionComparison, factorVariableInFormula, stepsize, nRepetitionsPerStepMax) {
    
    checkAndLoadPackages(c("tidyverse", "checkmate"), verbose = FALSE)
    
    # Check if contrasts have been specified correctly
    conditionsContrast = strsplit(conditionComparison, split = ",", fixed = TRUE)[[1]]
    assertVector(conditionsContrast, len = 2)
    
    assertFileExists(file_sampleTable, access = "r")
    
    assertIntegerish(stepsize, lower = 1)
    assertIntegerish(nRepetitionsPerStepMax, lower = 1)
    
    sampleData.df = read_tidyverse_wrapper(file_sampleTable, type = "tsv", col_names = TRUE) 
    
    assertSubset(conditionsContrast, sampleData.df$conditionSummary)
    
    # Get original ratio
    
    conditionCounter = table(sampleData.df$conditionSummary)
    
    nSamplesRareCondition     = min(conditionCounter)
    nSamplesFrequentCondition = max(conditionCounter)
    ratio = nSamplesRareCondition / (nSamplesFrequentCondition + nSamplesRareCondition)
    
    nameRareCondition      = names(conditionCounter)[conditionCounter == min(conditionCounter)]
    nameFrequentCondition  = names(conditionCounter)[conditionCounter == max(conditionCounter)]
    indexRareCondition     = which(sampleData.df$conditionSummary == nameRareCondition)
    indexFrequentCondition = which(sampleData.df$conditionSummary == nameFrequentCondition)
    
    stopifnot(length(indexRareCondition) == nSamplesRareCondition)
    
    if (nSamplesRareCondition > 10) {
        samplesRareCondition = c(2:5, seq(10, nSamplesRareCondition, stepsize))
        samplesRareCondition = c(3:9, seq(10, nSamplesRareCondition, stepsize))
        
        if (nSamplesRareCondition %% stepsize != 0) {
            samplesRareCondition = c(samplesRareCondition, nSamplesRareCondition)
        }
    } else {
        samplesRareCondition = c(2:nSamplesRareCondition)
    }
    
    nSamplesBase = length(samplesRareCondition)
    
    samplesFrequentCondition = ceiling(samplesRareCondition * (1/ratio - 1))
    
    # Correct rounding errors
    if (samplesFrequentCondition[nSamplesBase] > nSamplesFrequentCondition) samplesFrequentCondition[nSamplesBase] = nSamplesFrequentCondition
    
    # for each particular number of samples for the rare case
    
    subsamples.l = list()
    for (sampleBaseCur in 1:nSamplesBase) {
        
        nValidSamples = 0
        nPermutations = 0
        while (nValidSamples < nRepetitionsPerStepMax || nPermutations > 100) {
            
            nPermutations = nPermutations + 1
            # 1. Rare Condition 
            # How many different samples are actually possible?
            nCombinations = choose(nSamplesRareCondition, samplesRareCondition[sampleBaseCur])
            nSamplesCur = min(nCombinations, nRepetitionsPerStepMax)
            
            
            # Generate them
            table.l = list()
            sampleCombinations.l = list()
            for (i in 1:nSamplesCur) {
                
                table.l[[i]] = sample_n(sampleData.df[indexRareCondition,], samplesRareCondition[sampleBaseCur], replace = FALSE)
                
            }
            
            
            # 2. Frequent Condition 
            
            # How many different samples are actually possible?
            nCombinations = choose(nSamplesFrequentCondition, samplesFrequentCondition[sampleBaseCur])
            nSamplesCur = min(nCombinations, nRepetitionsPerStepMax)
            
            # Generate them
            table2.l = list()
            for (i in 1:nSamplesCur) {
                table2.l[[i]] = sample_n(sampleData.df[indexFrequentCondition,], samplesFrequentCondition[sampleBaseCur], replace = FALSE)
                
            }
            
            # Merge
            for (i in 1:min(length(table.l), length(table2.l))) {
                
                if (nValidSamples == nRepetitionsPerStepMax) {
                    break
                }
                
                listname = paste0(samplesRareCondition[sampleBaseCur], "_", i)
                subsamples.l[[listname]]  = bind_rows(table.l[[i]], table2.l[[i]])
                
                
                # Check validity of sample
                currentPermutationSampleComb = sort(unique(subsamples.l[[listname]]$SampleID))
                nUnique = length(unique(unlist(subsamples.l[[listname]][,factorVariableInFormula])))
                
                if (any(sapply(sampleCombinations.l, function(x) {identical(currentPermutationSampleComb, x)} )) || nUnique == 1) {
                    
                    # sample invalid, redo
                    sampleValid = FALSE
                } else {
                    
                    # sample unique and ok
                    sampleCombinations.l[[i]] = currentPermutationSampleComb
                    
                    sampleValid = TRUE
                    nValidSamples = nValidSamples + 1
                }
                
            } # end merge
            
        } # end while not enough valid samples
        
        
    }
    
    subsamples.l
}

testExistanceAndCreateDirectoriesRecursively <- function(directories) {
    
    for (dirname in unique(directories)) {
        
        if (!testDirectoryExists(dirname)) {
            dir.create(dirname, recursive = TRUE)
        } else {
            assertDirectoryExists(dirname, access = "w")
        }
        
    }
}

###########
# DESeq 2 #
###########

convertToFormula <- function(userFormula, validColnames = NULL) {
    
    # Create formula based on user-defined design
    designFormula = tryCatch({
        as.formula(userFormula)
    }, warning = function(w) {
        stop("Converting the design formula \"", userFormula, "\" created a warning, which should be checked carefully.")
    }, error = function(e) {
        stop("Design formula \"", userFormula, "\" not valid")
    })
    
    
    # Check colmn names
    if (!is.null(colnames)) {
        formulaVariables = attr(terms(designFormula), "term.labels")
        
        indexInteractionTerms = which(grepl(":", formulaVariables))
        if (length(indexInteractionTerms) > 0) {
            
            interactionTerms = formulaVariables[indexInteractionTerms]
            # Check whether all individual items have been specified
            for (interactionTermCur in interactionTerms) {
                componentsCur = strsplit(interactionTermCur, ":")[[1]]
                if (!all(componentsCur %in% formulaVariables)) {
                    message = paste0("Design formula is incorrect, not all terms from the interactions (", paste0(componentsCur, collapse = ","), ") have been specified individually.")
                    checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
                    
                }
                
            }
            formulaVariables = formulaVariables[-indexInteractionTerms]
            
            
        }
        
        assertSubset(formulaVariables, validColnames)
    }
    
    
    designFormula
    
}


glog2 <- function(x) ((asinh(x) - log(2))/log(2))

# Credits to Bernd Klaus
myMAPlot <- function(M, idx, main, minMean = 0) {
    M <- na.exclude(M[, idx])
    M <- (glog2(M))
    mean <- rowMeans(M)
    M <- M[mean > minMean, ]
    mean <- mean[mean > minMean]
    
    difference <- M[, 2] - M[, 1]
    
    main <- paste(main, ", Number of genes:", dim(M)[1]) 
    
    pl <- (qplot(mean, difference, main = main, ylim = c(-5,5), asp = 2, geom = "point", alpha = I(.5), color = I("grey30"), shape = I(16))  #  
           + geom_hline(aes(yintercept = 0), col = "#9850C3", show.legend = FALSE)
           + geom_smooth(method = "loess", se = FALSE, col = "#5D84C5", span = .4)
           + theme_bw()
    )
    
    return(pl)
}

.prettyNum <- function(number, verbose = TRUE) {
    prettyNum(number, big.mark = ",", scientific = FALSE)
}


.determineRowsVector <- function(nRows) {
    
    nRows = nrow(vsd)
    if (nRows < 500) {
        rowsVec = c(nRows)
    } else if (nRows < 5000) {
        rowsVec = c(500, nRows)
    } else if (nRows < 50000) {
        rowsVec = c(500, 5000, nRows)
    } else  {
        rowsVec = c(500, 5000, 50000, nRows)
    }
    
    rowsVec
}


splitStringInMultipleLines <- function(input, width, sepChar = "\n", verbose = TRUE) {
    
    assertCharacter(input, min.len = 1)
    assertInt(width, lower = 1)
    assertCharacter(sepChar, len = 1)
    
    
    as.character(sapply(input, function(x) {
        paste0(strsplit(x, paste0("(?<=.{", width, "})"), perl = TRUE)[[1]], collapse = sepChar)
    }
    ))
    
}

shortenStringsDataFrame <- function(df, width) {
    
    dataTypes = sapply(df, class)
    for (colCur in 1:ncol(df)) {
        
        if (dataTypes[colCur] %in% c("character", "factor")) {
            df[, colCur] = splitStringInMultipleLines(as.character(unlist(df[,colCur])), width = 20, sepChar = "\n", verbose = FALSE)
        }
    }
    
    df
}


plotPCAWrapper <- function(dd, PCAVariables, file = NULL, ...) {
    
    start = Sys.time()
    flog.info(paste0("Plotting PCA plots", if_else(is.null(file), "", paste0(" to file ", file))))
    
    checkAndLoadPackages(c("tidyverse", "checkmate", "DESeq2", "RColorBrewer"), verbose = FALSE)
    
    # Add newlines for columns with long values
    colData(dd) = shortenStringsDataFrame(colData(dd), 20)
   
    vsd <- DESeq2::varianceStabilizingTransformation(dd, blind=TRUE)

    metadata = colData(dd)
   
    rownames(metadata) = metadata$SampleID
    stopifnot(identical(colData(vsd)$SampleID, rownames(metadata)))
    
    
    assertSubset(PCAVariables, colnames(metadata))

    if (!is.null(file)) {
        pdf(file, ...)
    }
    
    rowsVec = .determineRowsVector(nrow(vsd))
    
    
    
    plots.l = list()
    
    for (nrowsCur in rowsVec) {
        
        index1 = paste0("top", nrowsCur)
        plots.l[[index1]] = list()
        
        for (varCur in PCAVariables) {
            
            index2 = paste0(varCur)
            
            pcadata = DESeq2::plotPCA(vsd, intgroup = varCur, ntop = nrowsCur, returnData = TRUE)
            
            if (nrowsCur == nrow(vsd)) {
                title = paste0("All ", .prettyNum(nrowsCur), " genes")
            } else {
                title = paste0("Top ", .prettyNum(nrowsCur), " genes")
            }
            
            
            nCategories = length(unique(unlist(metadata[,varCur])))
            if (nCategories > 8) {
                colors = colorRampPalette(RColorBrewer::brewer.pal(8, "Accent"))(nCategories)
            } else {
                colors = RColorBrewer::brewer.pal(8, "Accent")
            }
            
            percentVar <- round(100 * attr(pcadata, "percentVar"))
            g = ggplot(pcadata, aes_string("PC1", "PC2", color = varCur)) +
                geom_point(size = 4) + 
                xlab(paste0("PC1: ",percentVar[1],"% variance")) +
                ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
                coord_fixed() + theme_bw() + ggtitle(title) + theme(plot.title = element_text(hjust = 0.5)) +
                theme(legend.position = "bottom", legend.text = element_text(size = 7))
            
            if (is.factor(pcadata[,varCur])  || is.character(pcadata[,varCur])) {
                g = g + scale_color_manual(values = colors)
                
            } else {
                # let it draw colors automatically for now
            }
            
            plot(g)
            plots.l[[index1]][[index2]] = g
            
        }
        
    }
    
    if (!is.null(file)) {
        dev.off()
    }
    
    .printExecutionTime(start)
    
    plots.l
    
}

# PCA vs covariates adj. R-square plot
# from Mikael
plotPCACovariates <- function(df, metadata, variables, nPCs = 10, file, ...) {
    
    assertInteger(nPCs, lower = 1, upper = 10)
    assertSubset(variables, colnames(metadata))

    start = Sys.time()
    flog.info(paste0("PCA vs covariates adj. R-square", if_else(is.null(file), "", paste0(" to file ", file))))
    
    checkAndLoadPackages(c("tidyverse", "reshape2"), verbose = FALSE)
    
    
    if (!is.null(file)) {
        pdf(file, ...)
    }

    metadata = metadata[,variables]
    # Eliminate columns with less than 2 distinct values 
    dataTypes = sapply(metadata, class)
    deleteColumns = c()
    for (colCur in 1:ncol(metadata)) {
        
        if (dataTypes[colCur] %in% c("character", "factor")) {
            allValues = unique(as.character(unlist(metadata[, colCur])))
            naValues = which(allValues == "NA" | allValues == "na")
            if (length(naValues) > 0) allValues = allValues[-naValues]
            if (length(allValues) < 2) {
                flog.info(paste0("Ignore variable ", colnames(metadata)[colCur], " because it has fewer than 2 distinct non-NA values"))
                deleteColumns = c(deleteColumns, colCur)
            }
            
        }
    }
    metadata = metadata[, -deleteColumns]
    scale <- FALSE ##prcomp is used for PCA, should your data be scaled?

    
    r2_fun <- function(x, y) {
        summary(lm(y ~ x))$adj.r.squared
    } #function to get adj. R2 values
    
    rv <- rowVars(as.matrix(df))
    
    rowsVec = .determineRowsVector(nrow(df))
    for (nFeatures in rowsVec) {

        #select top X variables peaks
        select <- order(rv, decreasing = TRUE)[seq_len(nFeatures)]  
        pc <- prcomp(t(df[select,]), scale=scale)
        pcs_cv <- pc$x[,1:nPCs]   #Extract selected PCs
        
        eig.val <- 100*summary(pc)$importance[2,1:nPCs]
        
        labels <- paste0(colnames(pcs_cv), "\n",  " (", round(eig.val, 1),  " %)")
        
        stopifnot(identical(rownames(metadata), colnames(df)))

        
        r2_cv <- matrix(NA, nrow = ncol(metadata), ncol = ncol(pcs_cv),
                        dimnames = list(colnames(metadata), colnames(pcs_cv)))
 
        #get adj. R2 value for covariates ~ PCs
        for (cov in colnames(metadata)) {
            for (pc in colnames(pcs_cv)) {
                r2_cv[cov, pc] <- r2_fun(unlist(metadata[,cov]), pcs_cv[,pc])
            }
        }      
        
        # Remove NA only rows
        allNAs = apply(r2_cv, 1 , function(x) all(is.na(x)) )
        if (length(which(allNAs)) > 0) {
            flog.info(paste0("Ignore variables ", paste0(names(which(allNAs)), collapse = ","), " because they contained only non-finite values"))
            r2_cv = r2_cv[-which(allNAs),]
        }
        
        g = ggplot(reshape2::melt(r2_cv), aes(Var2, Var1, fill = value)) + geom_tile() +
            labs(x = "", y = "") + 
            scale_x_discrete(expand = c(0, 0), labels = labels) +
            scale_y_discrete(expand = c(0, 0)) +
            scale_fill_gradient2(limits = c(-1,1)) +
            labs(fill = "adj. R-sqr") +
            ggtitle(paste0("PCA for top ", nFeatures ," variable peaks"))
        print(g)
    }
    
    
    if (!is.null(file)) {
        dev.off()
    }
    
    .printExecutionTime(start)
    
}


plotHeatmaps <- function(dd, useBlind = TRUE, variables = NULL, file = NULL, ...) {
    
    start = Sys.time()
    flog.info(paste0("Plotting heatmaps", if_else(is.null(file), "", paste0(" to file ", file))))
    
    checkAndLoadPackages(c("tidyverse", "checkmate", "DESeq2", "RColorBrewer", "pheatmap"), verbose = FALSE)


    if (!is.null(file)) {
        pdf(file, ...)
    }
    
    #######################
    # Clustering heatmaps #
    #######################
    
    flog.info(paste(" Transform data, this may take a while ..."))
    
    # In order to test for differential expression, we operate on raw counts and use discrete distributions as described in the previous section on differential expression. However for other downstream analyses – e.g. for visualization or clustering – it might be useful to work with transformed versions of the count data. Maybe the most obvious choice of transformation is the logarithm. Since count values for a gene can be zero in some conditions (and non-zero in others), some advocate the use of pseudocounts, i.e. transformations of the form y=log2(n+n0) where n represents the count values and n0 is a positive constant.
    
    # In this section, we discuss two alternative approaches that offer more theoretical justification and a rational way of choosing parameters equivalent to n0 above. One makes use of the concept of variance stabilizing transformations (VST) (Tibshirani 1988; Huber et al. 2003; Anders and Huber 2010), and the other is the regularized logarithm or rlog, which incorporates a prior on the sample differences (Love, Huber, and Anders 2014). Both transformations produce transformed data on the log2 scale which has been normalized with respect to library size or other normalization factors.
    
    # The point of these two transformations, the VST and the rlog, is to remove the dependence of the variance on the mean, particularly the high variance of the logarithm of count data when the mean is low. Both VST and rlog use the experiment-wide trend of variance over mean, in order to transform the data to remove the experiment-wide trend. Note that we do not require or desire that all the genes have exactly the same variance after transformation. Indeed, in a figure below, you will see that after the transformations the genes with the same mean do not have exactly the same standard deviations, but that the experiment-wide trend has flattened. It is those genes with row variance above the trend which will allow us to cluster samples into interesting groups.
    
    
    
    ntd <- normTransform(dd)
    vsd <- vst(dd, blind=useBlind)
    rld <- rlog(dd, blind=useBlind)
    df <- as.data.frame(colData(dd))
    
    annotationCol = NULL
    if (!is.null(variables)) {
        df = df[,variables]
        annotationCol = df
    }
    
    
    dd.transform = varianceStabilizingTransformation(dd, blind = useBlind)
    
    
    flog.info(paste(" Plot data..."))
    
    rowsVec = .determineRowsVector(nrow(dd))
    
    for (nrowsCur in rowsVec) {
        
        select <- order(rowMeans(counts(dd,normalized=TRUE)),
                        decreasing=TRUE)[1:nrowsCur]
        
        ###################
        # 1. Count matrix #
        ###################
        
        mainLabel = paste0("After variance stabilizing transformation,\ntop ", nrowsCur, " genes w.r.t norm. counts")
        print(pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
                       cluster_cols=TRUE, annotation_col= annotationCol, main = mainLabel))
        
        mainLabel = paste0("After regularized log transformation\n, top ", nrowsCur, " genes w.r.t norm. counts")
        print(pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
                       cluster_cols=TRUE, annotation_col= annotationCol, main = mainLabel))
        
        mainLabel = paste0("After log2(n + 1) transformation\n, top ", nrowsCur, " genes w.r.t norm. counts")
        print(pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
                       cluster_cols=TRUE, annotation_col= annotationCol, main = mainLabel))
        
        
    }
    
    
    
    ######################################
    # 2. Heatmap of this distance matrix #
    ######################################
    # Another use of the transformed data is sample clustering. Here, we apply the dist function to the transpose of the transformed count matrix to get sample-to-sample distances.
    
    # A heatmap of this distance matrix gives us an overview over similarities and dissimilarities between samples. We have to provide a hierarchical clustering hc to the heatmap function based on the sample distances, or else the heatmap function would calculate a clustering based on the distances between the rows/columns of the distance matrix.
    
    sampleDists <- dist(t(assay(dd.transform)))
    
    sampleDistMatrix <- as.matrix(sampleDists)
    rownames(sampleDistMatrix) <- paste(dd.transform$SampleID, dd.transform$conditionSummary, sep="-")
    colnames(sampleDistMatrix) <- NULL
    colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
    
    ######################################
    # 2. Heatmap of this distance matrix #
    ######################################
    print(pheatmap(sampleDistMatrix,
                   clustering_distance_rows=sampleDists,
                   clustering_distance_cols=sampleDists,
                   col = colors, 
                   scale = "none",
                   main = "Based on sample distance matrix"))
    
    

    if (!is.null(file)) {
        dev.off()
    }
    
    .printExecutionTime(start)
    
    
}

createDebugFile <- function(snakemake) {
    
    checkAndLoadPackages(c("checkmate", "tools", "futile.logger"), verbose = FALSE)
    
    if (!testClass(snakemake, "Snakemake")) {
        flog.warn(paste0("Could not find snakemake object, therefore not saving anyting."))
        return()
    }
    
    logfile = snakemake@log[[1]][1]
    
    assertCharacter(logfile, any.missing = FALSE)
    
    filename = paste0(tools::file_path_sans_ext(logfile), ".rds")
    
    # flog.info(paste0("Saved Snakemake object for manually rerunning the R script to ", filename))
    
    saveRDS(snakemake, filename)
    
}

#' @importFrom BiocParallel multicoreWorkers MulticoreParam
#' @import checkmate 
.initBiocParallel <- function(nWorkers, verbose = FALSE) {
    
    checkAndLoadPackages(c("BiocParallel"), verbose = verbose)  
    assertInt(nWorkers, lower = 1)    
    assertInt(multicoreWorkers())
    
    if (nWorkers > multicoreWorkers()) {
        warning("Requested ", nWorkers, " CPUs, but only ", multicoreWorkers(), " are available and can be used.")
        nWorkers = multicoreWorkers()
    }
    
    MulticoreParam(workers = nWorkers, stop.on.error = TRUE)
    
}

#' @importFrom BiocParallel bplapply bpok
#' @import checkmate 
.execInParallelGen <- function(nCores, returnAsList = TRUE, listNames = NULL, iteration, abortIfErrorParallel = TRUE, verbose = TRUE, functionName, ...) {
    
    checkAndLoadPackages(c("BiocParallel"), verbose = verbose)  
    start.time  <-  Sys.time()
    
    assertInt(nCores, lower = 1)
    assertFlag(returnAsList)
    assertFunction(functionName)
    assertVector(iteration, any.missing = FALSE, min.len = 1)
    assert(checkNull(listNames), checkCharacter(listNames, len = length(iteration)))
    
    res.l = list()
    
    if (nCores > 1) {
        
        res.l = tryCatch( {
            bplapply(iteration, functionName, ..., BPPARAM = .initBiocParallel(nCores))
            
        }#, error = function(e) {
        #     warning("An error occured while executing the function with multiple CPUs. Trying again using only only one CPU...")
        #     lapply(iteration, functionName, ...)
        # }
        )
        
        failedTasks = which(!bpok(res.l))
        if (length(failedTasks) > 0) {
            warning("At least one task failed while executing in parallel, attempting to rerun those that failed: ",res.l[[failedTasks[1]]])
            if (abortIfErrorParallel) stop()
            
            res.l = tryCatch( {
                bplapply(iteration, functionName, ..., BPPARAM = .initBiocParallel(nCores), BPREDO = res.l)
                
            }, error = function(e) {
                warning("Once again, an error occured while executing the function with multiple CPUs. Trying again using only only one CPU...")
                if (abortIfErrorParallel) stop()
                lapply(iteration, functionName, ...)
            }
            )
        }
        
    } else {
        res.l = lapply(iteration, functionName, ...)
    }
    
    end.time  <-  Sys.time()
    
    if (nCores > multicoreWorkers()) {
        nCores = multicoreWorkers()
    }
    
    flog.info(paste0(" Finished execution using ",nCores," cores. TOTAL RUNNING TIME: ", round(end.time - start.time, 1), " ", units(end.time - start.time),"\n"))
    
    
    if (!returnAsList) {
        return(unlist(res.l))
    }
    
    if (!is.null(listNames)) {
        names(res.l) = listNames
    }
    
    res.l
    
}

chooseTFLabelSize <- function(nTF_label) {
  case_when( nTF_label < 20 ~ 8,
             nTF_label < 40 ~ 7,
             nTF_label < 60 ~ 6,
             nTF_label < 80 ~ 5,
             nTF_label < 100 ~ 4,
             TRUE ~ 3)
}

.printExecutionTime <- function(startTime) {
    
    endTime  <-  Sys.time()
    flog.info(paste0(" Finished sucessfully. Execution time: ", round(endTime - startTime, 1), " ", units(endTime - startTime)))
}

testFunction <- function(matrix) {
    
    matrix2 = matrix(data = 1, nrow = 10, ncol = 10)
    checkAndLogWarningsAndErrors(NULL, FALSE, isWarning = FALSE, dumpfile = paste0(getwd(), "/diffTF_error_dump.RData")) 
}


checkAndLogWarningsAndErrors <- function(object, checkResult, isWarning = FALSE, dumpfile = paste0(getwd(), "/diffTF_error_dump.RData")) {
    
    assert(checkCharacter(checkResult, len = 1), checkLogical(checkResult))
    
    if (checkResult != TRUE) {
        
        objectPart = ""
        if (!is.null(object)) {
            objectname = deparse(substitute(object))
            objectPart = paste0("Assertion on variable \"", objectname, "\" failed: ")
        } 
        
        lastPartError   = "# An error occurred. See details above. If you think this is a bug, please contact us. #\n"
        hashesStrError = paste0(paste0(rep("#", nchar(lastPartError) - 1), collapse = ""), "\n")
        messageError    = paste0(objectPart, checkResult, "\n\n", hashesStrError, lastPartError, hashesStrError)
        
        lastPartWarning = "# This warning may or may not be ignored. Carefully check the diagnostic files and downstream results. #\n"
        hashesStrWarning = paste0(paste0(rep("#", nchar(lastPartWarning) - 1), collapse = ""), "\n")
        messageWarning  = paste0(objectPart, checkResult, "\n\n", hashesStrWarning, lastPartWarning, hashesStrWarning)
        
        
        
        if (isWarning) {
            flog.warn(messageWarning)
            warning(messageWarning)
        } else {
            
            flog.info(paste0("Saving a dump file to ", dumpfile, ". You may use this file to trouble-shoot or send to others."))
            save.image(file = dumpfile)
            
            flog.error(messageError)
            stop(messageError)
        }
    }
}


# Which transformation of the y values to do?
transform_yValues <- function(values, addPseudoCount = TRUE, nPermutations, onlyForZero = TRUE) {
    
    # Should only happen with the permutation-based approach
    zeros = which(values == 0)
    if (length(zeros) > 0 & addPseudoCount) {
        values[zeros] = 1 / nPermutations
    }
    
    -log10(values)
}

filterLowlyExpressedGenes <- function(dd, comparisonMode, minMeanGroup, minMedianAll) {
  
  assertIntegerish(minMeanGroup)
  assertIntegerish(minMedianAll) 
    
  if (comparisonMode == "pairwise") {
      
    dd_counts = DESeq2::counts(dd, normalize = FALSE)
    
    samples_cond1 = colData(dd)$SampleID[which(colData(dd)$conditionSummary == levels(colData(dd)$conditionSummary)[1])]
    samples_cond2 = colData(dd)$SampleID[which(colData(dd)$conditionSummary == levels(colData(dd)$conditionSummary)[2])]
    
    # If in either of the two conditions counts are below a minimum OR the median is 0, we discard the gene and mark them as non-expressed
    idx <- (rowMeans(dd_counts[,samples_cond1]) > minMeanGroup | rowMeans(dd_counts[,samples_cond2]) > minMeanGroup) & rowMedians(dd_counts) > minMedianAll
    
  } else {
    
    # Here, we have to employ a different approach for when to filter genes given that we may have more than two conditions
    idx <- rowMedians(dd_counts) > minMedianAll
  }
  
  nFiltered = length(which(idx == FALSE))
  nKept     = length(which(idx == TRUE))
  
  if (nFiltered > 0) {
    flog.info(paste0("Filtered ", nFiltered, " genes from RNA-Seq table because of low counts. Remaining: ", nKept ))
    dd = dd[idx,]
  } 
  
  dd
}



transform_yValues_caption <- function() {
    "-log10"
}

my.median = function(x) median(x, na.rm = TRUE)
my.mean   = function(x) mean(x, na.rm = TRUE)


checkDesignIntegrity <- function(snakemake, par.l, sampleData.df) {
    
    
    par.l$designFormulaVariableTypes = snakemake@config$par_general$designVariableTypes
    checkAndLogWarningsAndErrors(par.l$designFormulaVariableTypes, checkCharacter(par.l$designFormulaVariableTypes, len = 1, min.chars = 3))
    par.l$designFormulaVariableTypes = gsub(" ", "", par.l$designFormulaVariableTypes)
    components = strsplit(par.l$designFormulaVariableTypes, ",")[[1]]
    
    par.l$conditionComparison  = snakemake@config$par_general$conditionComparison
    checkAndLogWarningsAndErrors(par.l$conditionComparison, checkCharacter(par.l$conditionComparison, len = 1))
    
    par.l$designFormula= snakemake@config$par_general$designContrast
    designFormula = as.formula(par.l$designFormula)
    
    formulaVariables = attr(terms(designFormula), "term.labels")
    
    indexInteractionTerms = which(grepl(":", formulaVariables))
    if (length(indexInteractionTerms) > 0) {
        
        interactionTerms = formulaVariables[indexInteractionTerms]
        # Check whether all individual items have been specified
        for (interactionTermCur in interactionTerms) {
            componentsCur = strsplit(interactionTermCur, ":")[[1]]
            if (!all(componentsCur %in% formulaVariables)) {
                message = paste0("Design formula is incorrect, not all terms from the interactions (", paste0(componentsCur, collapse = ","), ") have been specified individually.")
                checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
                
            }
            
        }
        formulaVariables = formulaVariables[-indexInteractionTerms]
        
        
    }
    
    checkAndLogWarningsAndErrors(components, checkVector(components, min.len = length(formulaVariables)))
    
    # Extract the variable that defines the contrast. Always the last element in the formula
    variableToPermute = formulaVariables[length(formulaVariables)]
    
    # Split further
    components2 = strsplit(components, ":")
    
    if (!all(sapply(components2,length) == 2)) {
        
        message = "The parameter \"designVariableTypes\" has not been specified correctly. It must contain all the variables that appear in the parameter \"designContrast\". See the documentation for details"
        checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
        
    }
    
    components3 = unlist(lapply(components2, "[[", 1))
    components3types = tolower(unlist(lapply(components2, "[[", 2)))
    names(components3types) = components3
    checkAndLogWarningsAndErrors(formulaVariables, checkSubset(formulaVariables, components3))
    checkAndLogWarningsAndErrors(components3types, checkSubset(components3types, c("factor", "integer", "numeric", "logical")))
    
    # Check the sample table. Distinguish between the two modes: Factor and integer for conditionSummary
    if (components3types["conditionSummary"] == "logical" | components3types["conditionSummary"] == "factor") {
        
        nDistValues = length(unique(sampleData.df$conditionSummary))
        if (nDistValues != 2) {
            message = paste0("The column 'conditionSummary' must contain exactly 2 different values, but ", nDistValues, " were found.") 
            checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
        }
        
        conditionsVec = strsplit(par.l$conditionComparison, ",")[[1]]
        if (!testSubset(as.character(unique(sampleData.df$conditionSummary)), conditionsVec)) {
            message = paste0("The specified elements for the parameter conditionComparison (", par.l$conditionComparison,   ") do not correspond to what is specified in the sample summary table (", paste0(unique(sampleData.df$conditionSummary), collapse = ","), "). All elements of conditionComparison must be present in the column conditionSummary.")
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
    
    
    list(types = components3types, variableToPermute = variableToPermute)
    
}



#####################
# AR CLASSIFICATION #
#####################

# TODO: peakID
createBindingMatrixFromFiles <- function(mapping, peakCounts, rootOutdir, extensionSize, comparisonType) {
  
  start = Sys.time()
  flog.info(paste0("Create binary binding matrix across all TF from existing files"))
  
  # Loop through all TFs and determine whether or not a peak has a binding site for this TF, 
  # resulting in a binary matrix
  TF.peakMatrix.l = list()
  for (TFCur in mapping$HOCOID) {
    
    mapping.subset.df = subset(mapping, HOCOID == TFCur)
    gene.sel = unique(mapping.subset.df$ENSEMBL)
    if (length(gene.sel) > 1) {
      message = paste0("Mapping for ", TFCur, " not unique, take only the first mapping (", gene.sel[1], 
                       ") and discard the others (", paste0(gene.sel[-1], collapse = ","), ")")
      checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
      
    }
    
    fileCur = paste0(rootOutdir, "/TF-SPECIFIC/", TFCur, "/extension", 
                     extensionSize, "/", comparisonType, TFCur,  ".output.tsv.gz")
    # TF.output.df = read.table(file = fileCur, header = TRUE)
    TF.output.df = read_tidyverse_wrapper(fileCur, type = "tsv", col_types = cols())
    
    TF.peakMatrix.l[[TFCur]] = peakCounts$peakID %in% TF.output.df$peakID
  }

  # This is the peak (rows) and TF binding sites (columns)
  TF.peakMatrix.df = as.data.frame(TF.peakMatrix.l)
  
  # Sanity check
  if (all(rowSums(TF.peakMatrix.df) == 0)) {
    message = paste0("All counts are 0, something is wrong.")
    checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
  }
  
  .printExecutionTime(start)
  TF.peakMatrix.df
  
}



readHOCOMOCOTable <- function(file, delim = " ") {

  HOCOMOCO_mapping.df = read_tidyverse_wrapper(file, type = "delim", delim = delim, col_types = cols()) %>%
    mutate(ENSEMBL = gsub("\\..+", "", ENSEMBL, perl = TRUE)) # Clean ENSEMBL IDs
  
  assertSubset(c("ENSEMBL", "HOCOID"), colnames(HOCOMOCO_mapping.df))
  
  
  HOCOMOCO_mapping.df
}



permute_rowsColumnWise <- function(input) {
  
  start = Sys.time()
  flog.info(paste0("Shuffle each column independently"))
  
  nRows = nrow(input)
  for (i in 1:ncol(input)) {
    input[,i] = (dplyr::select(input, i) %>% pull()) [sample(nRows)]
  }
  .printExecutionTime(start)
  input
}

filterPeaksByRowMeans <- function(peakCounts, TF.peakMatrix = NULL, minMean = 1, idColumn = "peakID") {
  
  start = Sys.time()
  flog.info(paste0("Filter peaks with a mean across samples of smaller than ", minMean))
  
  index_lastColumn = which(colnames(peakCounts) == idColumn)
  stopifnot(length(index_lastColumn) == 1)
  
  rowMeans2 = rowMeans(peakCounts[,-index_lastColumn])
  rowsToDelete = which(rowMeans2 < minMean)
  if (length(rowsToDelete) > 0) {
    flog.info(paste0("Removed ", length(rowsToDelete), " peaks out of ", nrow(peakCounts), 
                     " because they had a row mean smaller than ", minMean, "."))
    
    if (!is.null(TF.peakMatrix)) {
      # Filter these peaks also from the peakCount matrix  
      stopifnot(nrow(TF.peakMatrix) == nrow(peakCounts))
      TF.peakMatrix = TF.peakMatrix[-rowsToDelete,]
    }

    peakCounts = peakCounts[-rowsToDelete,]
  }
  
  .printExecutionTime(start)
  
  if (!is.null(TF.peakMatrix)) { 
      list(peakCounts = peakCounts, bindingMatrix = TF.peakMatrix)
  } else {
      peakCounts
  }
  
}




normalizeCounts <- function(rawCounts, method = "quantile", idColumn, removeCols = c(), returnDESeqObj = FALSE) {
    
    start = Sys.time()
    checkAndLoadPackages(c("tidyverse", "futile.logger", "checkmate", "preprocessCore"), verbose = FALSE)
    
    flog.info(paste0("Normalize counts. Method: ", method, ", ID column: ", idColumn))
    
    if (is.null(idColumn)) {
        rawCounts = as.data.frame(rawCounts) %>%
            rownames_to_column("ENSEMBL")
        idColumn = "ENSEMBL"
    } 
        
    ids = as.character(unlist(rawCounts[,idColumn]))
    
    # Test for additional character columns 
    colTypes = sapply(rawCounts, class)
    colTypes_rm = colTypes[which(colTypes == "character" | colTypes == "factor" | 
                                     names(colTypes) %in% removeCols)]
    
    if (length(colTypes_rm) > 1) {
        colnames_rm = names(colTypes_rm)[which(names(colTypes_rm) != idColumn)]
        flog.info(paste0("Remove the following columns because they do not represent counts: ", paste0(colnames_rm, collapse = ",")))
        colnames_rm = c(colnames_rm, idColumn)
    } else {
        colnames_rm = c(idColumn)
    }
    
    if (length(colnames_rm) > 0) {
        colnames_samples = colnames(rawCounts)[-which(colnames(rawCounts) %in% colnames_rm)]
    } else {
        colnames_samples = colnames(rawCounts)
    }
    
    rmCols = which(colnames(rawCounts) %in% colnames_rm)
    
    if (method == "quantile") {
        
        if (length(rmCols) > 0) {
            input = as.matrix(rawCounts[,-rmCols])
        } else {
            input = as.matrix(rawCounts)
        }
        counts.norm = normalize.quantiles(input)
        
    } else if (method == "DESeq_sizeFactor") {
        
        if (length(rmCols) > 0) {
            sampleData.df = data.frame( sampleID = colnames(rawCounts)[-rmCols])
            countDataNew = as.data.frame(rawCounts[, -rmCols])
        } else {
            sampleData.df = data.frame( sampleID = colnames(rawCounts))
            countDataNew = as.data.frame(rawCounts)
        }

        rownames(countDataNew) = ids
        
        stopifnot(identical(sampleData.df$sampleID, colnames(countDataNew)))
        
        dd <- DESeqDataSetFromMatrix(countData = countDataNew,
                                     colData = sampleData.df,
                                     design = as.formula(" ~ 1"))
        
        dd = estimateSizeFactors(dd)
        counts.norm = DESeq2::counts(dd, normalized = TRUE)
        
        if (returnDESeqObj) {
            return(dd)
        }
        
        
    } else if (method == "none") {
        
        if (length(rmCols) > 0) {
            counts.norm = rawCounts[,-rmCols]
        } else {
            counts.norm = rawCounts
        }

    } else  {
        stop("Not implemented yet")
    }
    
    .printExecutionTime(start)
    
    counts.norm %>% 
        as.data.frame()  %>% 
        as_tibble() %>% 
        dplyr::mutate(!!as.name(idColumn) := ids) %>%
        dplyr::select(!!as.name(idColumn), everything()) %>%
        purrr::set_names(c(idColumn, colnames_samples))
    
}




intersectData <- function(countsRNA, countsATAC, idColumn_RNA = "ENSEMBL", idColumn_ATAC = "peakID") {
    
    start = Sys.time()
    flog.info(paste0("Subset RNA and ATAC and keep only shared samples"))
    
    stopifnot(idColumn_RNA %in% colnames(countsRNA))
    stopifnot(idColumn_ATAC %in% colnames(countsATAC))
    
    flog.info(paste0(" Number of samples for RNA before filtering: " , ncol(countsRNA) - 1))
    flog.info(paste0(" Number of samples for ATAC before filtering: ", ncol(countsATAC) - 1))
    
    # Subset ATAC and RNA to the same set of samples
    sharedColumns = intersect(colnames(countsRNA), colnames(countsATAC))
    
    if (length(sharedColumns) == 0) {
        message = "RNA and ATAC counts have no shared samples. Verify that the colum names in the RNA-seq counts file are identical to the names in the sample table."
        checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
    } else {
        
        countsRNA.df  = dplyr::select(countsRNA,  c(idColumn_RNA,  sharedColumns))
        countsATAC.df = dplyr::select(countsATAC, c(idColumn_ATAC, sharedColumns))
        
        flog.info(paste0(" ", length(sharedColumns), " samples (", paste0(sharedColumns, collapse = ","), ") are shared between the ATAC-Seq and RNA-Seq data"))
        
        notIntersecting_ATAC = setdiff(colnames(countsATAC),c(sharedColumns,idColumn_ATAC))
        notIntersecting_RNA  = setdiff(colnames(countsRNA),c(sharedColumns,idColumn_RNA))
        if (length(notIntersecting_ATAC) > 0 ) {
            flog.warn(paste0("The following samples from the ATAC-Seq data will be ignored for the classification due to missing overlap with RNA-Seq: ", paste0(notIntersecting_ATAC, collapse = ",")))
        }
        if (length(notIntersecting_RNA) > 0) {
            flog.warn(paste0("The following samples from the RNA-Seq data will be ignored for the classification due to missing overlap with ATAC-Seq: ", paste0(notIntersecting_RNA, collapse = ",")))
        }
    }
    
    flog.info(paste0(" Number of samples for RNA after filtering: " , ncol(countsRNA.df) - 1))
    flog.info(paste0(" Number of samples for ATAC after filtering: ", ncol(countsATAC.df) - 1))
    .printExecutionTime(start)
    
    list(RNA = countsRNA.df, ATAC = countsATAC.df)
}

filterHOCOMOCOTable <- function(HOCOMOCO_table, TFs) {
    
    HOCOMOCO_mapping.df.overlap <- HOCOMOCO_table  %>%
        dplyr::filter(HOCOID %in% TFs) %>%
        dplyr::distinct(ENSEMBL, HOCOID)
    
    if (nrow(HOCOMOCO_mapping.df.overlap) == 0) {
        message = paste0("Number of rows of HOCOMOCO_mapping.df.overlap is 0. Something is wrong with the mapping table or the filtering")
        checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
    }
    
    HOCOMOCO_mapping.df.overlap
}

# Filter ATAC counts to consensus peaks
# Obselete
filterATACByConsensusPeaks <- function(countsATAC, consensusPeaks) {
    
    start = Sys.time()
    flog.info(paste0("Filter ATAC-Seq counts by consensus peaks"))
    flog.info(paste0(" Number of rows in consensus peaks: " , nrow(consensusPeaks)))
    flog.info(paste0(" Number of rows in ATAC-seq counts: " , nrow(countsATAC)))
    countsATAC.filt = dplyr::filter(countsATAC, peakID %in% consensusPeaks$peakID)
    flog.info(paste0(" Number of rows in ATAC-seq counts after filtering: " , nrow(countsATAC.filt)))
    
    diff = nrow(countsATAC) - nrow(countsATAC.filt)
    if (diff > 0) {
        
        message = paste0(" Filtered ", diff, " rows in ATAC-Seq counts because they do not overlap with consensus peaks. Note that the matching is done only by peakID, not by exact positions!")
        checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
        
        if (nrow(countsATAC.filt) == 0) {
            message = paste0(" No rows left in ATAC-Seq counts. Something is wrong with the IDs. Make sure IDs match between consensus peaks and ATAC-Seq counts")
            checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
        }
    }
    
    .printExecutionTime(start)
    countsATAC.filt
}




correlateATAC_RNA <- function(countsRNA, countsATAC, HOCOMOCO_mapping, corMethod = "pearson") {

    start = Sys.time()
    flog.info(paste0("Correlate RNA-Seq and ATAC-Seq counts"))
    
    # Filter to only the TFs
    # In addition, the no of TF because multiple TFs can map to the same gene/ ENSEMBL ID
    # Also filter 0 count genes because they otherwise cause errors downstream
    rowMeans = rowMeans(dplyr::select(countsRNA, -ENSEMBL))
    countsRNA.norm.TFs.df = dplyr::filter(countsRNA, ENSEMBL %in% HOCOMOCO_mapping$ENSEMBL, rowMeans > 0)
    
    diff = nrow(countsRNA) - nrow(countsRNA.norm.TFs.df)
    if (diff > 0) {
      flog.info(paste0(" Retain ", nrow(countsRNA.norm.TFs.df), " rows from RNA data (filter non-TF genes and TF genes with 0 counts throughout and keep only unique ENSEMBL IDs)."))
    }
    
    if (nrow(countsRNA.norm.TFs.df) == 0) {
        message = " No rows remaining from RNA_Seq counts after filtering against ENSEMBL IDs from HOCOMOCO. Check your ENSEMBL IDs for overlap with the HOCOMOCO translation table."
        checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
    }
    
    HOCOMOCO_mapping.exp = dplyr::filter(HOCOMOCO_mapping, ENSEMBL %in% countsRNA.norm.TFs.df$ENSEMBL)
    flog.info(paste0(" Correlate RNA-Seq and ATAC-Seq counts for ", nrow(countsATAC), " peaks and ", nrow(countsRNA.norm.TFs.df), " unique TF genes."))
    flog.info(paste0(" Note: For subsequent steps, the same gene may be associated with multiple TF, depending on the translation table."))
    # Correlate TF gene counts with ATAC-Seq counts 
    # countsRNA:  rows: all TF genes, columns: all samples
    # countsATAC: rows: peak IDs, columns: samples
    # Transverse both for the cor function then
    
    # counts for ATAC peaks may be 0 throughout, then a warning is thrown
    
    cor.m = t(cor(t(dplyr::select(countsRNA.norm.TFs.df, -ENSEMBL)), t(dplyr::select(countsATAC, -peakID)), method = corMethod))
    
    colnames(cor.m) = countsRNA.norm.TFs.df$ENSEMBL
    rownames(cor.m) = countsATAC$peakID
    
    # Some entries in the HOCOMOCO mapping can be repeated (i.e., the same ID for two different TFs, such as ZBTB4.S and ZBTB4.D)
    # Originally, we deleted these rows from the mapping and took the first entry only
    # However, since TFs with the same ENSEMBL ID can still be different with respect to their TFBS, we now duplicate such genes also in the correlation table
    #HOCOMOCO_mapping.exp = HOCOMOCO_mapping.exp[!duplicated(HOCOMOCO_mapping.exp[, c("ENSEMBL")]),]
    #assertSubset(as.character(HOCOMOCO_mapping.exp$ENSEMBL), colnames(sort.cor.m))
    
    sort.cor.m = cor.m[,names(sort(colMeans(cor.m)))] 
    # Change the column names from ENSEMBL ID to TF names. 
    # Reorder to make sure the order is the same. Due to the duplication ID issue, the number of columns may increase after the column selection
    sort.cor.m = sort.cor.m[,as.character(HOCOMOCO_mapping.exp$ENSEMBL)] 
    colnames(sort.cor.m) = as.character(HOCOMOCO_mapping.exp$HOCOID)

    .printExecutionTime(start)
    sort.cor.m
}

computeForegroundAndBackgroundMatrices <- function(peakMatrix, sort.cor.m) {
    
    start = Sys.time()
    flog.info(paste0("Compute foreground and background as well as their median values per TF"))
    # TODO: Extend with a few steps before even
    
    # This binary matrix has peaks as rows and TFs as columns of whether or not a particular peak has a TFBS from this TF or not
    # The unique prevents colnames from changeing, with ".1" being added to it automatically in case of duplicate column names
    sel.TF.peakMatrix.df = peakMatrix[,unique(colnames(sort.cor.m))]
    
    
    # 1. Focus on peaks with TFBS overlaps
    t.cor.sel.matrix = sort.cor.m
    # Transform 0 values into NA for the matrix to speed up subsequent analysis
    t.cor.sel.matrix[(sel.TF.peakMatrix.df == 0)] = NA
    # Goal: Eliminate all correlation values for cases in which a peak has no TFBS for the particular TF
    # Same dimensions as the two matrices used for input.
    # Matrix multiplication here essentially means we only multiply each individual entry with either 1 or NA
    # The result is a matrix full of NAs and the remaining entries are the correlation values for peaks that have a TFBS for the particular TF
    t.cor.sel.matrix = sel.TF.peakMatrix.df * t.cor.sel.matrix
    # Gives one value per TF, designating the median correlation per TF across all peaks
    median.cor.tfs = sort(apply(t.cor.sel.matrix, MARGIN = 2, FUN = my.median))
    
    # 2. Background
    # Start with the same correlation matrix
    t.cor.sel.matrix.non = sort.cor.m
    # Transform 1 values into NA for the matrix to speed up subsequent analysis
    t.cor.sel.matrix.non[(sel.TF.peakMatrix.df == 1)] = NA
    t.cor.sel.matrix.non = sel.TF.peakMatrix.df + t.cor.sel.matrix.non
    # Gives one value per TF, designating the median correlation per TF
    median.cor.tfs.non <- sort(apply(t.cor.sel.matrix.non, MARGIN=2, FUN = my.median))
    
    # Not used thereafter
    median.cor.tfs.rest <- sort(median.cor.tfs - median.cor.tfs.non[names(median.cor.tfs)])
    
    .printExecutionTime(start)
    
    list(median_foreground = median.cor.tfs, 
         median_background = median.cor.tfs.non[names(median.cor.tfs)],
         foreground = t.cor.sel.matrix,
         background = t.cor.sel.matrix.non
    )
}



calculate_classificationThresholds <- function(background, par.l) {
    
    start = Sys.time()
    flog.info(paste0("Calculate classification thresholds"))
    act.rep.thres.l = list()
    for (thresholdCur in par.l$allClassificationThresholds) {
        
        act.rep.cur = quantile(sort(apply(background, MARGIN = 2, FUN = my.median)), 
                               probs = c(thresholdCur, 1-thresholdCur))
        # Enforce the thresholds to be at least 0, so we never have an activator despite a negative median correlation 
        # and a repressor despite a positive one
        act.rep.thres.l[[as.character(thresholdCur)]][1] = min(0, act.rep.cur[1])
        act.rep.thres.l[[as.character(thresholdCur)]][2] = max(0, act.rep.cur[2])
        
        flog.info(paste0(" Thresholds for repressor/activator for threshold ", thresholdCur, ": ", 
                         act.rep.thres.l[[as.character(thresholdCur)]][1], 
                         " and ", 
                         act.rep.thres.l[[as.character(thresholdCur)]][2]))
        
    }
    
    .printExecutionTime(start)
    act.rep.thres.l
}


finalizeClassificationAndAppend <- function(output.global.TFs, median.cor.tfs, act.rep.thres.l, par.l, t.cor.sel.matrix, t.cor.sel.matrix.non, significanceThreshold_Wilcoxon = 0.05) {
    
   checkAndLoadPackages(c("progress"), verbose = FALSE)
  
    start = Sys.time()
    flog.info(paste0("Finalize classification"))
    colnameMedianCor         = paste0("median.cor.tfs")
    colnameClassificationPVal= paste0("classification_distr_rawP")
    
    AR.data = as.data.frame(median.cor.tfs)
    AR.data$TF = rownames(AR.data)
    colnames(AR.data)[1] = colnameMedianCor
    
    output.global.TFs[,colnameMedianCor] = NULL
    
    if ("TF" %in% colnames(output.global.TFs)) {
      output.global.TFs = merge(output.global.TFs, AR.data, by = "TF",all.x = TRUE)
    } else if ("TF.name" %in% colnames(output.global.TFs)) {
      output.global.TFs = merge(output.global.TFs, AR.data, by = "TF.name",all.x = TRUE)
    } else {
      message = paste0("Could npt find column for merging.")
      checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
    }
    
    
    
    # Define classes, tidyverse style
    for (thresCur in names(act.rep.thres.l)) {
        thresCur.v = act.rep.thres.l[[thresCur]]
        colnameClassificationCur = paste0("classification_q", thresCur)
        output.global.TFs[, colnameClassificationCur] = case_when( is.na(output.global.TFs[,colnameMedianCor]) ~ "not-expressed",
                                                                   output.global.TFs[,colnameMedianCor] < thresCur.v[1] ~ "repressor",
                                                                   output.global.TFs[,colnameMedianCor] > thresCur.v[2] ~ "activator",
                                                                   TRUE ~ "undetermined")
        output.global.TFs[, colnameClassificationCur] = factor(output.global.TFs[, colnameClassificationCur], levels = names(par.l$colorCategories))
        
    }
    
    if (!is.null(significanceThreshold_Wilcoxon)) {
        
        assertNumber(significanceThreshold_Wilcoxon, lower = 0, upper = 1)
        flog.info(paste0(" Perform Wilcoxon test for each TF. This may take a few minutes."))
        output.global.TFs[,colnameClassificationPVal] = NULL
        # Do a Wilcoxon test for each TF as a 2nd filtering criterion
        
        pb <- progress_bar$new(total = ncol(t.cor.sel.matrix))
   
        for (TFCur in colnames(t.cor.sel.matrix)) {
            
            pb$tick()
            rowNo = which(output.global.TFs$TF == TFCur)
            
            # Should normally not happen
            if (length(rowNo) != 1) {
                
                save(list = ls(),                 file = paste0(getwd(), "/diffTF_error_dump_local.RData") )
                save(list = ls(all.names = TRUE), file = paste0(getwd(), "/diffTF_error_dump_global.RData") )
                
                message = paste0("Mismatch detected between TF names in the correlation matrix and the output table. Error occured for the TF ", TFCur, ". This should not happen. Contact the authors.")
                checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
                
            }
            
            # Removing NAs actually makes a difference, as these are "artifical" anyway here due to the two matrices let's remove them
            dataMotif      = na.omit(t.cor.sel.matrix[,TFCur])
            dataBackground = na.omit(t.cor.sel.matrix.non[,TFCur])
            
            # Skip if NA for median correlation
            if (is.na(output.global.TFs[rowNo, colnameMedianCor])) {
                output.global.TFs[rowNo,colnameClassificationPVal] = NA
                next
            }
            
            # Test the distributions
            if (output.global.TFs[rowNo, colnameMedianCor] > 0) {
                alternativeTest = "greater"
            } else {
                alternativeTest = "less"
            }
            
            testResults = wilcox.test(dataMotif, dataBackground, alternative = alternativeTest)
            
            stopifnot(length(rowNo) == 1)
            output.global.TFs[rowNo,colnameClassificationPVal] = testResults$p.value
            
        }
        
        
        ################################################
        # POST-FILTER: CHANGE SOME TFs TO UNDETERMINED #
        ################################################
        
        # Change the classification with the p-value from the distribution test
        

        
        for (thresholdCur in par.l$allClassificationThresholds) {
            
            flog.info(paste0("  Post filter: Doing Wilcoxon test for threshold ", thresholdCur))
            
            colnameClassification      = paste0("classification_q", thresholdCur)
            colnameClassificationFinal = paste0("classification_q", thresholdCur, "_final")
            
            output.global.TFs[,colnameClassificationFinal] = output.global.TFs[,colnameClassification]
            
            TFs_to_change = dplyr::filter(output.global.TFs, (!!as.name(colnameClassification) == "activator" | 
                                                                  !!as.name(colnameClassification) == "repressor") & 
                                              !!as.name( colnameClassificationPVal) > !!significanceThreshold_Wilcoxon) %>%
                            pull(TF)
            
            # Filter some TFs to be undetermined
            if (length(TFs_to_change) > 0) {
                flog.info(paste0("  Changing the following TFs to 'undetermined' because they were classified as either activator or repressor before but the Wilcoxon test was not significant: ", paste0(TFs_to_change, collapse = ",")))
                
                output.global.TFs[which(output.global.TFs$TF  %in% TFs_to_change), colnameClassificationFinal] = "undetermined"
            }
            
            
        } # enhd for each threshold
        
    } # end if doWilcoxon
    
    # Print a summary of the classification
    flog.info(" Summary of classification:")
    colnamesIndex = which(grepl("final", colnames(output.global.TFs)))
    for (colnameCur in colnamesIndex) {
        flog.info(paste0("  Column ", colnames(output.global.TFs)[colnameCur]))
        tbl = table(output.global.TFs %>% pull(colnameCur))
        flog.info(paste0("   ", paste0(names(tbl), ": ", tbl), collapse = ", "))
    }
    
    .printExecutionTime(start)
    output.global.TFs
    
}


# Density plots for TFs

plot_density <- function(foreground.m, background.m, file = NULL, ...) {
    
    start = Sys.time()
    flog.info(paste0("Plotting density plots with foreground and background for each TF", if_else(is.null(file), "", paste0(" to file ", file))))
    stopifnot(identical(colnames(foreground.m), colnames(background.m)))
    
    # 1. Determine maximum y-values across all TFs
    yMax = 2
    for (colCur in seq_len(ncol(foreground.m))) {
      
       n_notNA = length(which(!is.na(foreground.m[,colCur])))
       if (n_notNA > 1){
         yMaxCur = max(c(density(foreground.m[,colCur], na.rm = TRUE)$y, density(background.m[,colCur], na.rm = T)$y))
         
         if (yMaxCur > yMax) {
           yMax = yMaxCur + 0.1
         }
       }
       
    }
    
    if (!is.null(file)) {
        pdf(file, ...)
    }
    
    for (colCur in seq_len(ncol(foreground.m))) {
        
        TFCur = colnames(foreground.m)[colCur]
        dataMotif      = foreground.m[,colCur]
        dataBackground = background.m[,colCur]
        mainLabel = paste0(TFCur," (#TFBS = ",length(which(!is.na(dataMotif)))," )")
        
        n_notNA1 = length(which(!is.na(dataMotif)))
        n_notNA2 = length(which(!is.na(dataBackground)))
        if (n_notNA1 > 1 & n_notNA2 > 1 ){
          
          plot(density(dataMotif, na.rm = TRUE), xlim = c(-1,1), ylim = c(0,yMax),
               main= mainLabel, lwd=2.5, col="red", axes = FALSE, xlab = "Pearson correlation")
          abline(v=0, col="black", lty=2)
          legend("topleft",box.col = adjustcolor("white",alpha.f = 0), legend = c("Motif","Non-motif"), lwd = c(2,2),cex = 0.8, col = c("red","darkgrey"), lty = c(1,1) )
          axis(side = 1, lwd = 1, line = 0)
          axis(side = 2, lwd = 1, line = 0, las = 1)
          
          lines(density(dataBackground, na.rm = T), lwd = 2.5, col = "darkgrey")
          
        } else {
          flog.warn(paste0(" Not enough data for estimating densities for TF ", TFCur, ", skip for plotting."))
        }
        
        
    } 
    
    if (!is.null(file)) {
        dev.off()
    }
    
    .printExecutionTime(start)
    
}

plot_classCorrelations <- function(DESeq_obj, output.global.TFs, HOCOMOCO_mapping, par.l, file = NULL) {
  
  start = Sys.time()
  flog.info(paste0("Plot class correlations", if_else(is.null(file), "", paste0(" to file ", file))))

  
  if (!is.null(file)) {
    pdf(file, ...)
  }

  DESeq_results  = DESeq2::results(DESeq_obj) %>% 
      as.data.frame() %>% 
      rownames_to_column("ENSEMBL") %>% 
      as_tibble()
  
  # Prepare the RNASeq data
  TF.specific = left_join(HOCOMOCO_mapping, DESeq_results, by = "ENSEMBL") %>% 
    dplyr::filter(!is.na(baseMean))
  
  if (nrow(TF.specific) == 0) {
    message = "The Ensembl IDs from the translation table do not match with the IDs from the RNA-seq counts table."
    checkAndLogWarningsAndErrors(NULL, message, isWarning = FALSE)
  }

  for (thresholdCur in par.l$allClassificationThresholds) {
    
    colnameClassificationCur = paste0("classification_q", thresholdCur, "_final")
    
    output.global.TFs.merged = output.global.TFs %>%
      dplyr::filter(!!as.name(colnameClassificationCur) != "not-expressed")  %>%
      full_join(TF.specific, by = c( "TF" = "HOCOID"))  %>%
      mutate(baseMeanNorm = (baseMean - min(baseMean, na.rm = TRUE)) / (max(baseMean, na.rm = TRUE) - min(baseMean, na.rm = TRUE)) + par.l$minPointSize)   %>%
      dplyr::filter(!is.na(!!as.name(colnameClassificationCur))) 
    
    
    for (classificationCur in unique(output.global.TFs.merged[,colnameClassificationCur])) {
      
      output.global.TFs.cur = dplyr::filter(output.global.TFs.merged, !!as.name(colnameClassificationCur) == classificationCur)
      
      cor.res.l = list()
      for (corMethodCur in c("pearson", "spearman")) {
        cor.res.l[[corMethodCur]] = cor.test(output.global.TFs.cur$weighted_meanDifference, 
                                             output.global.TFs.cur$log2FoldChange, 
                                             method = corMethodCur)
      }
      
      titleCur = paste0(classificationCur, ": R=", 
                        signif(cor.res.l[["pearson"]]$estimate, 2), "/", 
                        signif(cor.res.l[["spearman"]]$estimate, 2), ", p-value ", 
                        signif(cor.res.l[["pearson"]]$p.value,2),  "/", 
                        signif(cor.res.l[["spearman"]]$p.value,2), "\n(Pearson/Spearman, stringency: ", thresholdCur, ")")
      
      g = ggplot(output.global.TFs.cur, aes(weighted_meanDifference, log2FoldChange)) + geom_point(aes(size = baseMeanNorm)) + 
        geom_smooth(method = par.l$regressionMethod, color = par.l$colorCategories[classificationCur]) + 
        ggtitle(titleCur) + 
        ylab("log2 fold-change RNA-seq") + 
        theme_bw() + theme(plot.title = element_text(hjust = 0.5))
      plot(g)
      
    }
  }
  
  if (!is.null(file)) {
    dev.off()
  }
  
  .printExecutionTime(start)

}



plot_AR_thresholds  <- function(median.cor.tfs, median.cor.tfs.non, par.l, act.rep.thres.l, file = NULL, ...) {
    
    start = Sys.time()
    flog.info(paste0("Plotting AR summary plot", if_else(is.null(file), "", paste0(" to file ", file))))
    
    xlab="median pearson correlation (r)"
    ylab=""
    
    xlimMax = max(abs(c(range(median.cor.tfs), range(median.cor.tfs.non))))
    xlim =  c(-(xlimMax * 1.1), (xlimMax * 1.1))
    ylim = c(1,length(median.cor.tfs.non))
    
    nTF = length(median.cor.tfs)
    
    stopifnot(identical(names(median.cor.tfs), names(median.cor.tfs.non)))
    
    
    if (!is.null(file)) {
        pdf(file,  ...)
    }
    
    par(mfrow = c(1,1))
    
    for (thresCur in names(act.rep.thres.l)) {
        thresCur.v = act.rep.thres.l[[thresCur]]
        
        thresCur_upper = (1 - as.numeric(thresCur)) * 100
        thresCur_lower = as.numeric(thresCur) * 100
        
        mainCur = paste0("Stringency: ", thresCur)
        
        plot(median.cor.tfs.non, 1:nTF,
             xlim = xlim, ylim = ylim, main=mainCur, xlab=xlab, ylab=ylab,
             col = adjustcolor("darkgrey",alpha=1), pch = 16, cex = 0.5,axes = FALSE)
        points(median.cor.tfs, 1:nTF,
               pch=16,  cex = 0.5, 
               col= case_when(median.cor.tfs > thresCur.v[2] ~ par.l$colorCategories["activator"],
                              median.cor.tfs < thresCur.v[1] ~ par.l$colorCategories["repressor"],
                              TRUE ~ par.l$colorCategories["undetermined"])
        ) 
        

        
        text(x  = c((thresCur.v[1]),(thresCur.v[2])),
             y = c(nTF, nTF), pos = c(2,4),
             labels  = c(paste0(thresCur_lower, " percentile\n", "(", round(thresCur.v[1],5), ")"), 
                       paste0(thresCur_upper, " percentile\n", "(", round(thresCur.v[2],5), ")")),
             cex = 0.7, col = c("black","black"))
        abline(v = thresCur.v[1], col= par.l$colorCategories["repressor"])
        abline(v = thresCur.v[2], col= par.l$colorCategories["activator"])
        
        dataPoints = c(median.cor.tfs.non, median.cor.tfs)
        xAxisLimits = c(min(dataPoints) * 1.1, max(dataPoints) * 1.1)
        xAxisLimits[1] = max(-1, xAxisLimits[1])
        xAxisLimits[2] = min(1, xAxisLimits[2])
        
        # Set the limits dynamically
        if (xAxisLimits[1] > -0.1 | xAxisLimits[2] < 0.1) {
            defaultLimits = seq(-1,1,0.02)
        } else {
            defaultLimits = seq(-1,1,0.1)
        }
        
        defaultLimits = defaultLimits[-c(which(defaultLimits < 0 & defaultLimits < xAxisLimits[1] - 0.2), 
                                         which(defaultLimits > 0 & defaultLimits > xAxisLimits[2] + 0.2))]
        
  
        axis(side = 1, at = defaultLimits, lwd = 1, line = 0, cex = 1)
        
    }
    
    if (!is.null(file)) {
        dev.off()
    }
    
    .printExecutionTime(start)
    
}


# Code from Armando Reyes
plot_heatmapAR <- function(TF.peakMatrix.df, HOCOMOCO_mapping.df.exp, sort.cor.m, par.l, 
                           median.cor.tfs, median.cor.tfs.non, act.rep.thres.l, finalClassification = NULL,  file = NULL, ...) {
    
    start = Sys.time()
    flog.info(paste0("Plotting AR heatmap", if_else(is.null(file), "", paste0(" to file ", file))))
    
    
    missingGenes = which(!HOCOMOCO_mapping.df.exp$HOCOID %in% colnames(sort.cor.m))
    if (length(missingGenes) > 0) {
        HOCOMOCO_mapping.df.exp = dplyr::filter(HOCOMOCO_mapping.df.exp, HOCOID %in% colnames(sort.cor.m))
    }
    
    if (!is.null(file)) {
        pdf(file, ...)
    }
    
    cor.r.pearson.m <- sort.cor.m[,as.character(HOCOMOCO_mapping.df.exp$HOCOID)]
    
    stopifnot(identical(colnames(cor.r.pearson.m), as.character(HOCOMOCO_mapping.df.exp$HOCOID)))
    
    BREAKS = seq(-1,1,0.05)
    diffDensityMat = matrix(NA, nrow = ncol(cor.r.pearson.m), ncol = length(BREAKS) - 1)
    rownames(diffDensityMat) = HOCOMOCO_mapping.df.exp$HOCOID
    
    TF_Peak_all.m <- TF.peakMatrix.df
    TF_Peak.m <- TF_Peak_all.m
    
    for (i in 1:ncol(cor.r.pearson.m)) {
        TF = colnames(cor.r.pearson.m)[i]
        ## for the background, use all peaks
        h_noMotif = hist(cor.r.pearson.m[,TF][TF_Peak_all.m[,TF] == 0], breaks = BREAKS, plot = FALSE)
        ## for the foreground use only peaks with less than min_mot_n different TF motifs
        h_Motif   = hist(cor.r.pearson.m[,TF][TF_Peak.m[,TF]     != 0], breaks = BREAKS, plot = FALSE)
        diff_density = h_Motif$density - h_noMotif$density
        diffDensityMat[rownames(diffDensityMat) == TF, ] <- diff_density
    }
    diffDensityMat = diffDensityMat[!is.na(diffDensityMat[,1]),]
    colnames(diffDensityMat) = signif(h_Motif$mids,1)
    # quantile(diffDensityMat)
    
    ## check to what extent the number of TF motifs affects the density values
    n_min = if_else(colSums(TF_Peak.m) < nrow(TF_Peak.m),colSums(TF_Peak.m), nrow(TF_Peak.m) - colSums(TF_Peak.m))
    names(n_min) = HOCOMOCO_mapping.df.exp$HOCOID#[match(names(n_min), as.character(tf2ensg$ENSEMBL))]
    n_min <- sapply(split(n_min,names(n_min)),sum)
    
    # Make sure n_min and diffDenityMat are compatible because some NA rows may have been filtered out for diffDensityMat
    n_min <- n_min[rownames(diffDensityMat)]
    #quantile(n_min)
    remove_smallN = which(n_min < par.l$threshold_minNoTFBS_heatmap)
    cor(n_min[-remove_smallN], rowMax(diffDensityMat)[-remove_smallN], method = 'pearson')
    
    factorClassificationPlot <- sort(median.cor.tfs, decreasing = TRUE)
    diffDensityMat_Plot = diffDensityMat[match(names(factorClassificationPlot), rownames(diffDensityMat)), ]
    diffDensityMat_Plot = diffDensityMat_Plot[!is.na(rownames(diffDensityMat_Plot)),]
    annotation_rowDF = data.frame(median_diff = factorClassificationPlot[match(rownames(diffDensityMat_Plot), 
                                  names(factorClassificationPlot))])
    
    # Define the annotation row data frame with one column per threshold, with each TF colored according to its classification status
    anno_rowDF = data.frame(matrix(NA, nrow = nrow(diffDensityMat_Plot), ncol = 0))
    rownames(anno_rowDF) = rownames(diffDensityMat_Plot)
    annotation_colors = list()
    for (thresCur in names(act.rep.thres.l)) {

        nameCur = paste0(as.numeric(thresCur)*100, " / ", (1 - as.numeric(thresCur))*100, " %")
        colBreaks = unique(c((-1),
                             act.rep.thres.l[[thresCur]][1], 
                             act.rep.thres.l[[thresCur]][2],
                             1))

        anno_rowDF[,nameCur] = cut(annotation_rowDF$median_diff, breaks = colBreaks, labels = c("repressor", "undetermined", "activator"))
        
        if (is.null(finalClassification)) {
            # Plot original colors
            colors = c(par.l$colorCategories["repressor"],par.l$colorCategories["not-expressed"], par.l$colorCategories["activator"])
            names(colors) = levels(anno_rowDF[,nameCur])
            annotation_colors[[nameCur]] = colors
        } else {
            # Lighter colors because not final classification
            colors = c("#f18384", par.l$colorCategories["not-expressed"], "#aadaa8")
            names(colors) = levels(anno_rowDF[,nameCur])
            annotation_colors[[nameCur]] = colors
        }
        
        
        # Incorporate the provided final classification in here
        if (!is.null(finalClassification)) {
            
            colnameTable = paste0("classification_q", thresCur, "_final")
            
            assertDataFrame(finalClassification)
            assertSubset(c(colnameTable, "TF"), colnames(finalClassification))
            
            colnameAnno = paste0(nameCur, " final")
            matchTables = match(rownames(anno_rowDF), finalClassification$TF)
            anno_rowDF[,colnameAnno] = as.character(finalClassification[matchTables, colnameTable])
            
            # QC: Check if there are any impossible transitions
            nRows = nrow(anno_rowDF[which(anno_rowDF[,nameCur] == "undetermined" & anno_rowDF[,colnameAnno] == "activator"),])
            if (nRows > 0) {
                stop("Inconsistency deteced")
            }
            colors = c(par.l$colorCategories["repressor"],par.l$colorCategories["not-expressed"], par.l$colorCategories["activator"])
            
            names(colors) = levels(anno_rowDF[,nameCur])
            annotation_colors[[colnameAnno]] = colors
            
        }

    }

    labelMain = paste0(as.numeric(thresCur)*100, " / ", (1 - as.numeric(thresCur))*100, " %")
    labelMain = "Summary density heatmap for each TF\nand classifications across stringencies"
    pheatmap(diffDensityMat_Plot, cluster_rows = FALSE, cluster_cols = FALSE,
             fontsize_row = 1.25, scale = 'row' , fontsize_col = 10, fontsize = 8, labels_col = c(-1, -0.5, 0, 0.5, 1),
             annotation_row = anno_rowDF, annotation_colors = annotation_colors, annotation_legend = FALSE,  main = labelMain, legend = TRUE)

    
    if (!is.null(file)) {
        dev.off()
    }
    
    .printExecutionTime(start)
    
    
} # end function



plotDiagnosticPlots <- function(dd, conditionComparison, contrast = NULL, file = NULL, design = "conditionSummary", comparisonType = "pairwise", maxPairwiseComparisons = 5, plotMA = FALSE, alpha = 0.05, counts.raw = NULL, counts.norm = NULL) {
    
    #checkAndLoadPackages(c("tidyverse", "checkmate", "geneplotter", "DESeq2", "vsn", "RColorBrewer", "DEGreport", "ashr"), verbose = FALSE)
    checkAndLoadPackages(c("tidyverse", "checkmate", "geneplotter", "DESeq2", "vsn", "RColorBrewer", "ashr"), verbose = FALSE)
    flog.info(paste0("Plotting various diagnostic plots"))
    
    assert(testClass(dd, "MArrayLM"), testClass(dd, "DESeqDataSet"))
    
    assert(checkNull(file), checkDirectory(dirname(file), access = "w"))
    
    if (!is.null(file)) {
        pdf(file, width = 15, height = 10)
    }
    
    if (testClass(dd, "DESeqDataSet")) {
        
        res_DESeq = tryCatch( {
            results(dd)
            
        }, error = function(e) {
            FALSE
        }
        )
        
        if (!is.null(contrast)) {
            if (is.logical(res_DESeq)) {
                res = DESeq2::results(dd, contrast = contrast)
            } else {
                res = res_DESeq
            }
            coefCur = paste0(contrast[1], "_", contrast[2], "_vs_", contrast[3])
            assertSubset(coefCur, resultsNames(dd))
        } else {
            if (is.logical(res_DESeq)) {
                res = DESeq2::results(dd)
            } else {
                res = res_DESeq
            }
            coefCur = length(resultsNames(dd))
        }
        
        counts.raw  = DESeq2::counts(dd, normalized = FALSE)
        counts.norm = DESeq2::counts(dd, normalized = TRUE)
        
        design.df <- as.data.frame(colData(dd))
        
        nSamples = nrow(colData(dd))
        colors = colorRampPalette(brewer.pal(9, "Set1"))(nSamples)
        
    } else {
        
        nSamples = ncol(counts.raw)
        colors = colorRampPalette(brewer.pal(9, "Set1"))(nSamples)
        
    }
    
    
    if (plotMA) {
        
        if (testClass(dd, "MArrayLM")) {
            checkAndLoadPackages(c("limma"))
            assertVector(conditionComparison, len = 2)
            title = paste0("limma results\n", conditionComparison[1], " vs. ", conditionComparison[2])
            isSign = if_else(p.adjust(dd$p.value[,ncol(dd$p.value)], method = "BH") < alpha, paste0("sign. (BH, ", alpha, ")"), "not-significant")
            
            flog.info(paste0(" Plotting MA plots from limma..."))
            
            limma::plotMA(dd, main = title, status = isSign)
            
        } else {
            
            # DESeq2 specific MA plots
            flog.info(paste0(" Plotting MA plot from DESeq (1)..."))
            
            log2fcRanges = range(res$log2FoldChange)
            if (log2fcRanges[1] < -8) {
                flog.info(paste0(" y-limits for log2fc ranges capped at -8, some genes have higher values"))
                log2fcRanges[1] = -8
            }
            if (log2fcRanges[2] > 8) {
                flog.info(paste0(" y-limits for log2fc ranges capped at 8, some genes have higher values"))
                log2fcRanges[2] = 8
            }
            # 1. Regular MA plot:
            # shows the log2 fold changes attributable to a given variable over the mean of normalized counts for all the samples in the DESeqDataSet
            DESeq2::plotMA(res, main = paste0("Regular MA plot\n(", coefCur, ")"), ylim = log2fcRanges)
            
            # 2. MA plot based on shrunken log2 fold changes
            # Shrinkage of effect size (LFC estimates) is useful for visualization and ranking of genes. To shrink the LFC, we pass the dds object to the function  lfcShrink. Below we specify to use the apeglm method for effect size shrinkage (Zhu, Ibrahim, and Love 2018), which improves on the previous estimator.
            
            # It is more useful visualize the MA-plot for the shrunken log2 fold changes, which remove the noise associated with log2 fold changes from low count genes without requiring arbitrary filtering thresholds.
            flog.info(paste0(" Plotting MA plot from DESeq (2)..."))
            
            dd_LFC <- lfcShrink(dd, res = res, type = "ashr")
            
            log2fcRanges = range(dd_LFC$log2FoldChange)
            if (log2fcRanges[1] < -8) {
                flog.info(paste0(" y-limits for log2fc ranges capped at -8, some genes have higher values"))
                log2fcRanges[1] = -8
            }
            if (log2fcRanges[2] > 8) {
                flog.info(paste0(" y-limits for log2fc ranges capped at 8, some genes have higher values"))
                log2fcRanges[2] = 8
            }
            
            DESeq2::plotMA(dd_LFC, main = paste0("MA plot based on shrunken log2 fold changes\n(", coefCur, ")"), ylim = log2fcRanges)
            
        }
    }
    
    # 2. Densities of counts for the different samples. 
    # Since most of the genes are (heavily) affected by the experimental conditions, a succesful normalization will lead to overlapping densities
    
    xlabCur = "Mean log counts"
    flog.info(paste0(" Plotting multidensity plot..."))
    
    if (nSamples < 10) {
        
        multidensity(log(counts.raw + 0.5), xlab = xlabCur, main = "Non-normalized log counts", col = colors)
        multidensity(log(counts.norm + 0.5) , xlab = xlabCur, main = "Normalized log counts", col = colors) 
        
        multiecdf(log(counts.raw + 0.5), xlab = xlabCur, main = "Non-normalized log counts", col = colors)
        multiecdf(log(counts.norm + 0.5) , xlab = xlabCur, main = "Normalized log counts", col = colors) 
        
    } else {
        
        warning("Omitting legend due to large number of samples (threshold is 10 at the moment)")
        
        multidensity(log(counts.raw + 0.5),  xlab = xlabCur, main = "Non-normalized log counts", legend = NULL, col = colors)
        multidensity(log(counts.norm + 0.5), xlab = xlabCur, main = "Normalized log counts",     legend = NULL, col = colors)
        
        multiecdf(log(counts.raw + 0.5),   xlab = xlabCur, main = "Non-normalized log counts", legend = NULL, col = colors)
        multiecdf(log(counts.norm + 0.5) , xlab = xlabCur, main = "Normalized log counts",     legend = NULL, col = colors)
    }
    
    if (testClass(dd, "DESeqDataSet")) {
        
        if (maxPairwiseComparisons > 0) {
            
            flog.info(paste0(" Plotting pairwise sample comparisons. This may take a while."))
            
            # 3. Pairwise sample comparisons.
            # To further assess systematic differences between the samples, we can also plot pairwise mean–average plots: We plot the average of the log–transformed counts vs the fold change per gene for each of the sample pairs.
            MA.idx = t(combn(seq_len(dim(colData(dd))[1]), 2))
            
            if (nrow(MA.idx) > maxPairwiseComparisons) {
                
                flog.info(paste0("The number of pairwise comparisons to plot exceeds the current maximum of ", maxPairwiseComparisons, ". Only ", maxPairwiseComparisons, " pairwise comparisons will be shown in the PDF."))
                MA.idx.filt = MA.idx[1:maxPairwiseComparisons,, drop = FALSE]
                
            } else {
                MA.idx.filt = MA.idx
            }
            
            for (i in seq_along(MA.idx.filt[,1])) { 
                
                flog.info(paste0(" Plotting pairwise comparison ", i, " out of ", nrow(MA.idx.filt)))
                label = paste0(colnames(dd)[MA.idx.filt[i,1]], " vs ", colnames(dd)[MA.idx.filt[i,2]])
                suppressWarnings(print(myMAPlot(counts.norm, c(MA.idx[i,1], MA.idx.filt[i,2]), main =  label)))
            }
            
            # Show an empty page with a warning if plots have been omitted
            if (nrow(MA.idx) > nrow(MA.idx.filt)) {
                
                plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
                message = paste0("All remaining pairwise comparisons plots\nbetween samples have been omitted\nfor time and memory reasons.\nThe current maximum is set to ", maxPairwiseComparisons, ".")
                text(x = 0.5, y = 0.5, message, cex = 1.6, col = "red")
            }
            
        }
        
        # 4. Mean SD plot: Plot row standard deviations versus row means
        notAllZeroPeaks <- (rowSums(DESeq2::counts(dd)) > 0)
        
        suppressWarnings(meanSdPlot(assay(dd[notAllZeroPeaks,])))
        
        # if (comparisonType == "pairwise") {
        #     # Mean-Variance QC plots: p-value distribution gives an idea on how well you model is capturing the input data and as well whether it could be some problem for some set of genes. In general, you expect to have a flat distribution with peaks at 0 and 1. In this case, we add the mean count information to check if any set of genes are enriched in any specific p-value range.
        #     # Variation (dispersion) and average expression relationship shouldn’t be a factor among the differentially expressed genes. When plotting average mean and standard deviation, significant genes should be randomly distributed.
        #     g = degQC(counts.norm, design.df[[design]], pvalue = res[["pvalue"]])
        #     plot(g)
        # }
        # 
        # 
        # # Covariates effect on count data: Another important analysis to do if you have covariates is to calculate the correlation between PCs from PCA analysis to different variables you may think are affecting the gene expression. This is a toy example of how the function works with raw data, where clearly library size correlates with some of the PCs.
        # 
        # resCov <- degCovariates(log2(counts.norm + 0.5), colData(dd))
        # 
        # # Covariates correlation with metrics: Also, the correlation among covariates and metrics from the analysis can be tested. This is useful when the study has multiple variables, like in clinical trials. The following code will return a correlation table, and plot the correlation heatmap for all the covariates and metrics in a table.
        # 
        # cor <- degCorCov(colData(dd))
        # # 
        # # library(pheatmap)
        # # pheatmap(log2(counts.norm + 0.5)[1:10,], 
        # #          annotation_col = as.data.frame(colData(dd))[,1:4],
        # #          annotation_colors = degColors(colData(dd)[1:4],
        # #                                        con_values = c("white","red")
        # #          )
        # # )
        # 
        
    }
    
    
    if (!is.null(file)) {
        dev.off()
    }
}

setDebugMode <- function(debugMode) {
    
    #assertSubset(debugMode, c("TRUE", "FALSE", "True", "False", "true", "false", "", TRUE, FALSE), empty.ok = TRUE)
    if (is.null(debugMode)) {
        return(FALSE)
    } else {
        if (as.logical(debugMode) == TRUE) {
            flog.info(paste0("diffTF debug mode is enabled. R session files will be stored for this script. Use them to troubleshoot errors you get or send them to the authors for investigation upon being asked to do so."))
            return(TRUE)
        } else {
            return(FALSE)
        }
    }
   
}




