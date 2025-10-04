rm(list=ls())
closeAllConnections() 
graphics.off()
################################################################################
# As long as this following section is user-verified to be ready, press 'Source' button or select all and run using CMD+R/CTRL+R to run all lines in R
################################################################################
# Set the following options based on your need.
# Refer to the README.md for detailed documentation on the sampling, QC filtering, and donor-wise consensus functions.
sampling <- F                          # Set to T to enable random sampling of sequences

maxSampleSize <- 100                   # Maximum number of sequences/donors to sample from (if sampling is TRUE)

filterToggleList <- c(T,F)             # Vector to cycle through QC filter on/off options. 
# Set to c(T) if you only want analysis on filtered data. 
# Set to c(F) if you only want analysis on unfiltered data.

donorwiseConsensusToggleList = c(T,F)  # Vector to cycle through donor-wise consensus on/off options. 
# Set to c(T) if you only want analysis on donor-wise consensus data. 
# Set to c(F) if you only want analysis on individual sequence data.

# Get working directory
MyWD <- getwd()
setwd(MyWD)

# HXB2 reference genome (included with the .R script file in the same directory) 
HXB2GenomePath = paste0(MyWD,"/HXB2.fasta")

############################## Dependencies ####################################
# The following version of HIVprimertestR was built using R version "4.4.2"

# Check below listed package versions (install if required). 
# BiocManager >= "1.30.25"  # install.packages("BiocManager") #BiocManager is required to install Biostrings, msa, and pwalign.
# Biostrings  >= "2.74.1"   # BiocManager::install("Biostrings")
# dplyr       >= "1.1.4"    # install.packages("dplyr")
# msa         >= "1.38.0"   # BiocManager::install("msa")
# pwalign     >= "1.2.0"    # BiocManager::install("pwalign")
# tools       >= "4.4.2"    # No need to install- available with base R

packages <- c("Biostrings", "dplyr", "msa", "pwalign", "tools")

for (pkg in packages) {
  library(pkg, character.only = TRUE)
  cat(pkg, "version:", as.character(packageVersion(pkg)), "\n")
}

# The following sections are automated and do not require user edits unless users wish to adjust specific settings

############################## Utilities #######################################

# Sample function (If there are more than "maxSize" sequences, randomly sample down to "maxSize")
randomSample <- function(stringSet, maxSize = 100) {
  if (length(stringSet) <= maxSize) return(stringSet)
  sampledIndices <- sample(seq_along(stringSet), maxSize)
  return(stringSet[sampledIndices])
}

# Degap function
degap <- function(stringSet){
  stringSetDegapped <- DNAStringSet(gsub("-","",stringSet))
  stringSetDegapped <- DNAStringSet(gsub(" ","",stringSetDegapped))
  return(stringSetDegapped)
}

# Convert formatted .csv to .fasta file
csvToFasta <- function(dfs) {
  stringSets <- list()
  
  # Check for column count
  badNames <- names(dfs)[which(sapply(dfs, ncol) != 4)]
  if(is.null(badNames)){
    if (length(badNames) > 0) {
      stop(paste0("The following INPUT dataframes do not have exactly 4 columns: ", 
                  paste(badNames, collapse = ", "), 
                  "\nPlease maintain INPUT columns in the following format: GroupName\tDonorID\tSeqID\tSequence"))
    }
  }
  
  
  # Process each dataframe
  for (i in seq_along(dfs)) {
    df <- dfs[[i]]
    primerName <- names(dfs)[i]
    
    stringSet <- degap(DNAStringSet(df[,4]))
    
    names(stringSet) <- apply(df[, 1:3], 1, function(row) {
      cleaned <- gsub("\\.", "-", row)
      paste(cleaned, collapse = ".")
    })
    
    stringSets[[primerName]] <- stringSet
  }
  
  return(stringSets)
}

# Split a .fasta/.csv based on population
splitByPopulation <- function(stringSets, inputDir) {
  splitStringSets <- list()
  
  for (i in seq_along(stringSets)) {
    stringSet <- stringSets[[i]]
    primerName <- names(stringSets)[i]
    seqNames <- names(stringSet)
    
    popNames <- sapply(strsplit(seqNames, "\\."), `[`, 1)
    stringSetDegapped <- degap(stringSet)
    
    for (population in unique(popNames)) {
      populationStringSetSubset <- stringSetDegapped[popNames == population]
      stringSetName <- paste0(population, ".", primerName)
      outPath <- file.path(inputDir, paste0(stringSetName, ".fasta"))
      
      writeXStringSet(populationStringSetSubset, filepath = outPath)
      splitStringSets[[stringSetName]] <- populationStringSetSubset
    }
  }
  return(splitStringSets)
}

# Function to delete gap-only columns
removeGapOnlyColumns <- function(stringSet) {
  # Convert to matrix
  mat <- as.matrix(stringSet)
  
  # Keep columns that contain at least one non-gap character
  non_gap_columns <- apply(mat, 2, function(col) any(col != "-"))
  
  # Subset matrix and rebuild sequences
  filtered_mat <- mat[, non_gap_columns, drop = FALSE]
  filtered_seqs <- apply(filtered_mat, 1, paste0, collapse = "")
  
  # Return cleaned DNAStringSet with names preserved
  result <- DNAStringSet(filtered_seqs)
  names(result) <- names(stringSet)
  return(result)
}

# Sequence sort function
sortBySequenceNames <- function(msaObject) {
  msaObject <- DNAStringSet(msaObject)
  # Extract the sequence names
  seq_names <- names(msaObject)
  
  # Order the sequence names
  ordered_indices <- order(seq_names)
  
  # Return the sorted DNAStringSet
  return(msaObject[ordered_indices])
}

# Sort sequences by percent representation
sortByRepresentationPercent <- function(msaObject, top.10=F) {
  msaObject <- DNAStringSet(msaObject)  # Convert input to DNAStringSet
  
  # Extract the percent representation (It is the 2nd to last item after strsplit)
  percentRepresentation <- as.numeric(gsub("%","",sapply(names(msaObject), function(x) strsplit(x, "_")[[1]][length(strsplit(x, "_")[[1]])-1])))
  
  # Handle NA values safely by ordering with `na.last = TRUE`
  ordered_indices <- order(percentRepresentation, decreasing = TRUE, na.last = TRUE)
  
  # Get top 10 if requested, but only if enough sequences exist
  if (top.10 && length(ordered_indices) >= 10) {
    ordered_indices <- ordered_indices[1:10]
  }
  
  # Return the sorted DNAStringSet
  return(msaObject[ordered_indices])
}

# Split msa string set based on split point & reference sequence
splitMSA <- function(msaObject, splitPoint){
  return(c(fiveEnd=subseq(msaObject,start=1,end=splitPoint),threeEnd=subseq(msaObject,start=splitPoint+1,end=width(msaObject))))
}

# Remove empty sequences and filter out short/long sequences
cleanStringSet <- function(stringSet, reference, outputDir, filterToggle){
  refLength <- nchar(reference)
  stringSetDG <- degap(stringSet)
  
  # Remove empty sequences (i.e., sequences of length 0)
  emptyIndices <- which(stringSetDG=="")
  
  if(!length(emptyIndices)==0){
    stringSet.flt0 <- stringSetDG[-emptyIndices]
  } else {
    stringSet.flt0 <- stringSetDG
  }
  
  # Filter out short or long sequences only if filterToggle is "True"
  if(filterToggle){
    # 1. Remove short sequences
    if(length(which(as.integer(sapply(stringSetDG,nchar))<(0.5*refLength)))>=1){
      stringSet.flt1 <- stringSetDG[-which(nchar(stringSetDG)<(0.5*refLength))]
    } else {
      stringSet.flt1 <- stringSetDG
    }
    
    # 2. Remove long sequences
    if(length(which(as.integer(sapply(stringSetDG,nchar))>(1.5*refLength)))>=1){
      stringSet.flt2 <- stringSet.flt1[-which(nchar(stringSet.flt1)>(1.5*refLength))]
    } else {
      stringSet.flt2 <- stringSet.flt1
    }
  } else {
    stringSet.flt2 <- stringSet.flt0
  }
  
  # Return filtered & edited stringSet
  return(stringSet.flt2)
}

# Generate donor0-wise consensus
generateConsensus <- function(stringSet.RAW, threshold=0.25) {
  stringSet <- degap(stringSet.RAW)
  
  # Vectorized donor name extraction (Donor name assumed to be in 2nd column of INPUT file)
  donorList <- vapply(strsplit(names(stringSet), "\\."), function(parts) {
    if (length(parts) >= 2) parts[2] else NA_character_
  }, character(1))
  
  # Split sequences by donor
  donorMap <- split(stringSet, donorList)
  donorMap <- donorMap[!is.na(names(donorMap))]  # Remove NA donors
  
  # Generate consensus for each donor
  consensusList <- unlist(lapply(donorMap, function(subset) {
    if (length(subset) > 1) {
      alignment <- msa(degap(subset), method = "Muscle", type = "dna")
      paste0(consensusString(alignment,threshold=threshold)[[1]])
    } else {
      paste0(subset[[1]])  # Already a DNAStringSet of length 1
    }}))
  
  # Combine into DNAStringSet
  consensusStringSet <- DNAStringSet(consensusList)
  
  return(consensusStringSet)
}

# Combine 5' and 3' LTR regions
# Prefer 3' LTR over 5' and keep one sequence per donor per primer
combine5and3 <- function(Five,Three){
  Three.complement.Names <- setdiff(names(Five),names(Three))
  Combined <- Three
  Combined <- append(Combined,Five[which(names(Five)%in%Three.complement.Names)])
  return(Combined)
}

# Collapse sequences
collapseIdenticalSequences <- function(stringSet){
  degappedStringSet <- degap(stringSet)
  table <- table(degappedStringSet)
  
  headers <- paste0("n",table,".",names(table))
  collapsedStringSet <- DNAStringSet(names(table))
  names(collapsedStringSet) <- headers
  
  return(collapsedStringSet)
}

# Hamming distance function accounting for degenerate DNA bases
calculateHammingDistance <- function(reference, stringSet) {
  # Degenerate base dictionary
  degeneracyDict <- list(
    W = c("A", "T"), S = c("C", "G"), M = c("A", "C"), K = c("G", "T"),
    R = c("A", "G"), Y = c("C", "T"), B = c("C", "G", "T"), D = c("A", "G", "T"),
    H = c("A", "C", "T"), V = c("A", "C", "G"), N = c("A", "C", "G", "T")
  )
  
  # Expand the reference once
  ref_chars <- strsplit(as.character(reference), "")[[1]]
  n <- length(ref_chars)
  
  # Convert degenerate reference bases once
  ref_sets <- lapply(ref_chars, function(char) {
    if (char %in% names(degeneracyDict)) degeneracyDict[[char]] else char
  })
  
  # Output vector
  hamming_distances <- numeric(length(stringSet))
  names(hamming_distances) <- names(stringSet)
  
  for (i in seq_along(stringSet)) {
    seq_chars <- strsplit(as.character(stringSet[[i]]), "")[[1]]
    
    dist <- 0
    for (j in seq_len(n)) {
      ref_base <- ref_chars[j]
      seq_base <- seq_chars[j]
      
      # Skip if exact match
      if (seq_base == ref_base) next
      
      seq_set <- if (seq_base %in% names(degeneracyDict)) degeneracyDict[[seq_base]] else seq_base
      ref_set <- ref_sets[[j]]
      
      # Compute fractional mismatch
      if (is.character(ref_set) && is.character(seq_set)) {
        if (any(seq_set %in% ref_set)) {
          dist <- dist + (1 - length(intersect(seq_set, ref_set)) / length(seq_set))
        } else {
          dist <- dist + 1
        }
      }
    }
    
    hamming_distances[i] <- dist
  }
  
  return(data.frame(HammingDistanceWithDegeneracy = hamming_distances))
}

# Calculate hamming distances and return summary files
calculateHammingDistances <- function(stringSets, outputDir, population, primerSubset, filterToggle){
  # Initialize empty dataframe
  hammingDistResults <- data.frame(
    Target_Sequence = character(),
    Sample_Size = integer(),
    Primer_Name = character(),
    Primer_Target_Sequence = character(),
    Representation_Percent = numeric(),
    Full_HD = numeric(),
    FivePrimeHalf_HD = numeric(),
    ThreePrimeHalf_HD = numeric(),
    Terminal2nt_HD = numeric(),
    Definition.1_Qualified = logical(),
    Definition.2_Qualified = logical(),
    Definition.3_Qualified = logical(),
    stringsAsFactors = FALSE
  )
  
  # Use setNames to preserve mapping between primer and its matches
  map <- setNames(lapply(names(primerSubset), function(pn) grep(pn, names(stringSets))), names(primerSubset))
  
  # For each primer, calculate HD against associated viral sub-populations
  for(primerName in unique(names(primerSubset))){
    stringSet <- degap(stringSets[[unlist(map[[primerName]])]])
    # To keep the primer on the top in sorted msa objects, add 00_ in front of the primerName
    primer <- primerSubset[primerName]
    names(primer) <- paste0("00_",primerName)
    msaObject <- DNAStringSet(msa(c(primer,stringSet),method="Muscle",type = "dna"))
    msaObjectSorted <- sortBySequenceNames(msaObject)
    
    # To identify where gaps are in the msa
    reference <- paste0(msaObjectSorted[[1]])
    refLength <- nchar(gsub("-","",reference))
    splitPositions <- which(strsplit(reference, "")[[1]] != "-")
    
    # Split msa object into 5' and 3' halves based on number of non-gap characters
    splitStringSetsMiddle <- splitMSA(msaObject=msaObjectSorted, splitPoint=splitPositions[floor(refLength/2)])
    
    if(grepl("reverse",primerName,ignore.case = T)){
      # Flip 5' and 3' for reverse primers
      fivePrimeHalf <- splitStringSetsMiddle$threeEnd
      threePrimeHalf <- splitStringSetsMiddle$fiveEnd
      terminal2nts <- splitMSA(msaObject=msaObjectSorted, splitPoint=splitPositions[2])$fiveEnd
    } else { # If the primer name doesn't specify reverse, treat it as forward or probe
      fivePrimeHalf <- splitStringSetsMiddle$fiveEnd
      threePrimeHalf <- splitStringSetsMiddle$threeEnd
      terminal2nts <- splitMSA(msaObject=msaObjectSorted, splitPoint=splitPositions[floor(refLength-1)])$threeEnd
    }
    
    # Get substitution matrices
    subMatrix <- calculateHammingDistance(reference=msaObjectSorted[[1]],stringSet=msaObjectSorted[-1])
    fivePrimeSubMatrix <- calculateHammingDistance(reference=fivePrimeHalf[[1]],stringSet=fivePrimeHalf[-1])
    threePrimeSubMatrix <- calculateHammingDistance(reference=threePrimeHalf[[1]],stringSet=threePrimeHalf[-1])
    terminal2SubMatrix <- calculateHammingDistance(reference=terminal2nts[[1]],stringSet=terminal2nts[-1])
    
    if(filterToggle){
      # Identify sequences for which the Hamming Distance is greater than 50% of the primer length
      badSequenceIndices <- which(subMatrix[,1]>refLength*0.5)
    } else {
      badSequenceIndices <- c()
    }
    
    # Re-conduct MSA after trimming bad sequences without nonsense sequences 
    if(length(badSequenceIndices)!=0){
      
      for(index in badSequenceIndices){
        badSequence <- msaObjectSorted[index]
        name <- names(msaObjectSorted[index])
        alignment <- pairwiseAlignment(degap(reference), degap(paste0(badSequence)), type="overlap", substitutionMatrix=NULL, fuzzyMatrix=NULL,  gapOpening=10, gapExtension=4, scoreOnly=FALSE)
        if(paste0(subject(alignment))!=""){
          msaObjectSorted[index] <- subject(alignment)
          names(msaObjectSorted)[index] <- name
        }
      }
      
      msaObjectFiltered <- msa(degap(msaObjectSorted),method="Muscle",type = "dna")
      msaObjectFilteredSorted <- sortBySequenceNames(msaObjectFiltered)
      
      reference <- paste0(msaObjectFilteredSorted[[1]])
      refLength <- nchar(gsub("-","",reference))
      splitPositions <- which(strsplit(reference, "")[[1]] != "-")
      
      # Split msa object into first and second half
      splitFilteredStringSetsMiddle <- splitMSA(msaObject=msaObjectFilteredSorted, splitPoint=splitPositions[floor(refLength/2)])
      
      if(grepl("reverse",primerName,ignore.case = T)){
        # Flip 5' and 3' for reverse primers
        fivePrimeHalfFiltered <- splitFilteredStringSetsMiddle$threeEnd
        threePrimeHalfFiltered <- splitFilteredStringSetsMiddle$fiveEnd
        terminal2ntsFiltered <- splitMSA(msaObject=msaObjectFilteredSorted, splitPoint=splitPositions[2])$fiveEnd
      } else { # If the primer name doesn't specify reverse, treat it as forward or probe
        fivePrimeHalfFiltered <- splitFilteredStringSetsMiddle$fiveEnd
        threePrimeHalfFiltered <- splitFilteredStringSetsMiddle$threeEnd
        terminal2ntsFiltered <- splitMSA(msaObject=msaObjectFilteredSorted, splitPoint=splitPositions[floor(refLength-1)])$threeEnd
      }
      
      # Get substitution matrices
      subMatrix.V02 <- calculateHammingDistance(reference=msaObjectFilteredSorted[[1]],stringSet=msaObjectFilteredSorted[-1])
      fivePrimeSubMatrix.V02 <- calculateHammingDistance(reference=fivePrimeHalfFiltered[[1]],stringSet=fivePrimeHalfFiltered[-1])
      threePrimeSubMatrix.V02 <- calculateHammingDistance(reference=threePrimeHalfFiltered[[1]],stringSet=threePrimeHalfFiltered[-1])
      terminal2SubMatrix.V02 <- calculateHammingDistance(reference=terminal2ntsFiltered[[1]],stringSet=terminal2ntsFiltered[-1])
    } else {
      msaObjectFilteredSorted <- msaObjectSorted
      subMatrix.V02 <- subMatrix
      fivePrimeSubMatrix.V02 <- fivePrimeSubMatrix
      threePrimeSubMatrix.V02 <- threePrimeSubMatrix
      terminal2SubMatrix.V02 <- terminal2SubMatrix
    }
    
    # Change qualification status based on reference sequence type
    if(grepl("probe", primerName, ignore.case = T)){
      qualificationStatus <- subMatrix.V02[,1]<=1
      qualificationStatus.without3Mismatch <- qualificationStatus
    } else {
      qualificationStatus <- (fivePrimeSubMatrix.V02[,1]<=3 & threePrimeSubMatrix.V02[,1]<=1)
      qualificationStatus.without3Mismatch <- ((fivePrimeSubMatrix.V02[,1]<=3 & threePrimeSubMatrix.V02[,1]<=1) & terminal2SubMatrix.V02[,1]<=0)
    }
    
    sampleSizes <- as.numeric(sub("^n(\\d+)\\..*", "\\1",row.names(subMatrix.V02)))
    collapsedSeq <- sub("^n\\d+\\.", "",row.names(subMatrix.V02))
    identityStatus <- subMatrix.V02[,1] == 0
    representationPercent <- sampleSizes/sum(sampleSizes) * 100
    
    tempResults <- data.frame(
      Target_Sequence = collapsedSeq,
      Sample_Size = sampleSizes,
      Primer_Name = primerName,
      Primer_Target_Sequence = paste0(primer),
      Representation_Percent = representationPercent,
      Full_HD = subMatrix.V02[,1],
      FivePrimeHalf_HD = fivePrimeSubMatrix.V02[,1],
      ThreePrimeHalf_HD = threePrimeSubMatrix.V02[,1],
      Terminal2nt_HD = terminal2SubMatrix.V02[,1],
      Definition.1_Qualified = qualificationStatus,
      Definition.2_Qualified = qualificationStatus.without3Mismatch,
      Definition.3_Qualified = identityStatus,
      stringsAsFactors = FALSE
    )
    
    hammingDistResults <- rbind(hammingDistResults,tempResults)
    
    #Change names of sequences and export the fasta
    names(msaObjectFilteredSorted)[-1] <- paste0(population, "_", gsub("00_","",primerName), "_", 
                                                 round(tempResults$Representation_Percent, digits=2), "%_", 
                                                 row.names(subMatrix.V02))
    
    msaObjectSortedByRepresentation <- sortByRepresentationPercent(msaObjectFilteredSorted[-1])
    outputMSA <- c(msaObjectFilteredSorted[1],msaObjectSortedByRepresentation)
    writeXStringSet(outputMSA, filepath = paste0(outputDir,"Representation/",population,".",gsub("00_","",primerName),".MSA.fasta"))
  }
  
  hammingDistResults$population <- population
  
  write.csv(hammingDistResults,paste0(outputDir,"HammingDistanceRAW/",population,"_hammingDistance.csv"),row.names = F)
  
  HD.Summary.row.names <- c("SampleSize",
                            "Median HD(1st & 3rd Quantile) Full Length",
                            "Median HD(1st & 3rd Quantile)  5' Half",
                            "Median HD(1st & 3rd Quantile) 3' Half",
                            "Median HD(1st & 3rd Quantile) 2nts 3' terminal",
                            "PercentQualified.Definition.1", "PercentQualified.Definition.2", "PercentQualified.Definition.3")
  
  HDSummary <- data.frame(matrix(nrow=length(HD.Summary.row.names),ncol=length(primerSubset)))
  
  row.names(HDSummary)  <- HD.Summary.row.names
  colnames(HDSummary) <- gsub("00_","",names(primerSubset))
  
  for (primerName in names(primerSubset)) {
    subset <- hammingDistResults[gsub("00_","",hammingDistResults$Primer_Name) == gsub("00_","",primerName),]
    
    fullLength_HDs <- na.omit(rep(subset$Full_HD,subset$Sample_Size))
    fivePrime_HDs <- na.omit(rep(subset$FivePrimeHalf_HD,subset$Sample_Size))
    threePrime_HDs <- na.omit(rep(subset$ThreePrimeHalf_HD,subset$Sample_Size))
    terminal2nt_HDs <- na.omit(rep(subset$Terminal2nt_HD,subset$Sample_Size))
    
    fullLength_SummaryStats <- summary(fullLength_HDs)
    fivePrime_SummaryStats <- summary(fivePrime_HDs)
    threePrime_SummaryStats <- summary(threePrime_HDs)
    terminal2nt_SummaryStats <- summary(terminal2nt_HDs)
    
    HDSummary[[gsub("00_","",primerName)]] <- c(sum(subset$Sample_Size),
                                                paste0(fullLength_SummaryStats['Median']," (",fullLength_SummaryStats['1st Qu.']," to ",fullLength_SummaryStats['3rd Qu.'],")"),
                                                paste0(fivePrime_SummaryStats['Median']," (",fivePrime_SummaryStats['1st Qu.']," to ",fivePrime_SummaryStats['3rd Qu.'],")"),
                                                paste0(threePrime_SummaryStats['Median']," (",threePrime_SummaryStats['1st Qu.']," to ",threePrime_SummaryStats['3rd Qu.'],")"),
                                                paste0(terminal2nt_SummaryStats['Median']," (",terminal2nt_SummaryStats['1st Qu.']," to ",terminal2nt_SummaryStats['3rd Qu.'],")"),
                                                paste0(round(sum(subset[subset$Definition.1_Qualified,]$Representation_Percent),digits=2),"%"),
                                                paste0(round(sum(subset[subset$Definition.2_Qualified,]$Representation_Percent),digits=2),"%"),
                                                paste0(round(sum(subset[subset$Definition.3_Qualified,]$Representation_Percent),digits=2),"%"))
  }
  
  write.csv(HDSummary,paste0(outputDir,population,"_HDSummary.csv"))
  
  return(HDSummary)
}

# Function to run the pipeline
runFullAnalysis <- function(
    sampling = F,
    maxSampleSize = 100,
    filterToggleList = c(T,F),
    donorwiseConsensusToggleList = c(T,F),
    HXB2GenomePath = "./HXB2.fasta"
) {
  # Read HXB2 reference genome
  HXB2Genome <- readDNAStringSet(HXB2GenomePath)
  
  # Check if INPUT directory exists
  inputDir <- paste0(getwd(), "/INPUT/")
  if (!dir.exists(inputDir)) stop("INPUT directory not found!")
  
  # Step 1: Read input files
  
  # Step 1.1: Read PrimersINPUT.csv file
  primerPath <- paste0(inputDir,"PrimersINPUT.csv")
  if (!file.exists(primerPath)) {
    stop("Error: PrimersINPUT.csv file does not exist in the INPUT folder.")
  } else {
    primerDF <- read.csv(primerPath)
    primers <- DNAStringSet(primerDF[,2])
    names(primers) <- primerDF[,1]
  }
  
  # Step 1.2: Read & handle trimmed sequence query files
  InputFilepaths <- list.files(path = inputDir, pattern = "\\.csv$", full.names = TRUE, recursive = FALSE)
  queryTargetSetFilepaths <- InputFilepaths[!grepl("PrimersINPUT\\.csv",InputFilepaths)]
  
  if (!is.null(queryTargetSetFilepaths)) {
    if (length(queryTargetSetFilepaths) > 0) {
      csvDfs <- lapply(queryTargetSetFilepaths, read.csv)
      names(csvDfs) <- gsub("\\.csv","",basename(queryTargetSetFilepaths))
      queryTargetSets <- csvToFasta(csvDfs)
      queryTargetSetsSplitbyPopulation <- splitByPopulation(queryTargetSets,inputDir)
    }
  } else {
    stop("Error: No INPUT query datasets found!")
  }
  
  print("Input files prepared")
  
  # Step 2: (Optional) Apply sampling
  if (sampling && length(queryTargetSetsSplitbyPopulation) > 0) {
    queryTargetSetsSplitbyPopulation <- lapply(queryTargetSetsSplitbyPopulation, function(stringSet) {
      randomSample(stringSet = stringSet, maxSize = maxSampleSize)
    })
  }
  
  # Step 3: Apply filters and consensus logic
  resultsDir <- paste0(getwd(),"/RESULTS/")
  dir.create(resultsDir)
  
  for(filterToggle in filterToggleList){
    for(donorwiseConsensusToggle in donorwiseConsensusToggleList){
      print(paste0("Now running analysis with filter toggle ", filterToggle, " & donor-wise consensus toggle ", donorwiseConsensusToggle))
      
      # Create output directories
      outputDir <- paste0(resultsDir,"Filter.",filterToggle,"_DonorConsensus.",donorwiseConsensusToggle,"/")
      dirPaths <- file.path(outputDir, c("","Consensus", "Combined", "Collapsed", "Representation","HammingDistanceRAW"))
      sapply(dirPaths, dir.create, recursive = TRUE, showWarnings = FALSE)
      
      # Initiate empty lists to store results
      targetRegions.Filtered <- list()
      targetRegions.donorwiseConsensus <- list()
      
      # Step 1. Filter and remove empty sequences, fix names, and generate consensus
      for (stringSetName in names(queryTargetSetsSplitbyPopulation)) {
        stringSet <- queryTargetSetsSplitbyPopulation[[stringSetName]]
        splitName <- strsplit(stringSetName, "\\.")[[1]]
        primerName <- splitName[2]
        
        # Step 1: Filter & clean the datasets
        cleanedSet <- cleanStringSet(stringSet = stringSet,
                                     reference = primers[[primerName]],
                                     filterToggle = filterToggle)
        
        targetRegions.Filtered[[stringSetName]] <- cleanedSet
        
        # Step 2: Generate donor-wise consensus if toggled
        if (donorwiseConsensusToggle) {
          consensusSet <- generateConsensus(stringSet.RAW=cleanedSet,threshold=0.25)
          targetRegions.donorwiseConsensus[[stringSetName]] <- consensusSet
          
          outPath <- file.path(outputDir, "Consensus",
                               paste0(stringSetName, ".Consensus.fasta"))
          writeXStringSet(consensusSet, filepath = outPath)
        } else {
          targetRegions.donorwiseConsensus[[stringSetName]] <- cleanedSet
        }
      }
      
      print("Input sequences cleaned and donor-wise consensus sequences generated (if toggled)")
      
      # Step 3. Subset for each population 
      # Parse metadata
      file_info <- data.frame(
        stringSetName = names(targetRegions.donorwiseConsensus),
        stringsAsFactors = FALSE
      )
      
      file_info$population <- sapply(strsplit(file_info$stringSetName, "\\."), `[`, 1)
      
      populations <- unique(file_info$population)
      
      for (population in populations) {
        targetRegions.ConsensusCombined <- list()
        targetRegions.Collapsed <- list()
        
        popSubset <- file_info[file_info$population == population, ]
        
        print(paste0("For primer/probe: ",primerName,". Now generating results for group ", population))
        
        # Extract subset of DNAStringSets
        popNames <- popSubset$stringSetName
        populationSubset <- targetRegions.donorwiseConsensus[popNames]
        
        # Combine 5' and 3' LTR if present
        fiveLTRs <- popNames[grep("5LTR", popNames)]
        threeLTRs <- popNames[grep("3LTR", popNames)]
        ltrRoots <- unique(gsub("\\.(5|3)LTR", "", c(fiveLTRs, threeLTRs)))
        
        for (root in ltrRoots) {
          Five <- populationSubset[[paste0(root,".5LTR")]]
          Three <- populationSubset[[paste0(root,".3LTR")]]
          if (!is.null(Five) && !is.null(Three)) {
            targetRegions.ConsensusCombined[[root]] <- combine5and3(Five, Three)
          }
        }
        
        # Add non-LTR regions
        nonLTRNames <- setdiff(popNames, c(fiveLTRs, threeLTRs))
        for (name in nonLTRNames) {
          targetRegions.ConsensusCombined[[name]] <- populationSubset[[name]]
        }
        
        # Step 3: Write combined sequences and collapse
        for (name in names(targetRegions.ConsensusCombined)) {
          combinedPath <- file.path(outputDir, "Combined", paste0(name, ".Combined.fasta"))
          writeXStringSet(targetRegions.ConsensusCombined[[name]], filepath = combinedPath)
          
          collapsed <- collapseIdenticalSequences(targetRegions.ConsensusCombined[[name]])
          targetRegions.Collapsed[[name]] <- collapsed
          
          collapsedPath <- file.path(outputDir, "Collapsed", paste0(name, ".Collapsed.fasta"))
          writeXStringSet(collapsed, filepath = collapsedPath)
        }
        
        # Get primer names
        primerNames <- unique(sapply(strsplit(file_info$stringSetName, "\\."), `[`, 2))
        
        # Step 4: Calculate Hamming Distances
        HDSummary <- calculateHammingDistances(
          stringSets = targetRegions.Collapsed,
          outputDir = outputDir,
          population = population,
          primerSubset = primers[primerNames],
          filterToggle = filterToggle
        )
      }
    }
  }
  
  print("Analysis Finished!")
}

############################## Run pipeline ####################################

runFullAnalysis(sampling = sampling,
                maxSampleSize = maxSampleSize,
                filterToggleList = filterToggleList,
                donorwiseConsensusToggleList = donorwiseConsensusToggleList,
                HXB2GenomePath = HXB2GenomePath)
