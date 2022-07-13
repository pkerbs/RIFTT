addPSToTable <- function(rt){
  pl <- list()
  pl <- data.frame("gene"=unique(c(rt$gene1,rt$gene2)),
                         ar_num_partners=0,
                         fc_num_partners=0,stringsAsFactors = F)
  
  for(j in 1:nrow(pl)){
    partner <- c()
    temp <- AR_results
    partner <- c(partner,subset(temp,gene1==pl$gene[j])$gene2)
    partner <- c(partner,subset(temp,gene2==pl$gene[j])$gene1)
    partner <- unique(partner)
    pl$ar_num_partners[j] <- length(partner)
    partner <- c()
    temp <- FC_results
    partner <- c(partner,subset(temp,gene1==pl$gene[j])$gene2)
    partner <- c(partner,subset(temp,gene2==pl$gene[j])$gene1)
    partner <- unique(partner)
    pl$fc_num_partners[j] <- length(partner)
  }
  pl$PS <- (pl$fc_num_partners + pl$ar_num_partners)/2
    
  rt$PS <- -1
  for(i in 1:nrow(rt)){
    rt$PS[i] <- (pl$PS[pl$gene==rt$gene1[i]] + pl$PS[pl$gene==rt$gene2[i]])/2
  }
  return(rt)
}

addRSToTable <- function(rt,FTSpassIdx){
  # Compile recurrence library of raw calls
  RSlib <- AR_results[,match(c("sample","gene1","gene2"),colnames(AR_results))]
  RSlib$gene1 <- gsub("\\(\\d+\\)","",RSlib$gene1)
  RSlib$gene2 <- gsub("\\(\\d+\\)","",RSlib$gene2)
  RSlib$caller <- "AR"
  RSlib <- RSlib %>%
    group_by(gene1,gene2,caller) %>%
    summarize(rec=length(unique(sample)),.groups='drop')
  temp <- FC_results[,match(c("sample","gene1","gene2"),colnames(FC_results))]
  temp$caller <- "FC"
  temp <- temp %>%
    group_by(gene1,gene2,caller) %>%
    summarize(rec=length(unique(sample)),.groups='drop')
  RSlib <- rbind(RSlib,temp)
  
  # Compile recurrence library of calls passing the FTS filter
    FTSpassed <- resultTable[FTSpassIdx,match(c("sample","gene1","gene2","caller"),colnames(resultTable))]
    FTSpassed$gene1 <- gsub("\\(\\d+\\)","",FTSpassed$gene1)
    FTSpassed$gene2 <- gsub("\\(\\d+\\)","",FTSpassed$gene2)
    FTSpassed <- FTSpassed %>%
      group_by(gene1,gene2,caller) %>%
      summarize(rec=length(unique(sample)),.groups='drop')
  
  # Calculate the Robustness Score  
    rt$RS <- 0
    for(i in 1:nrow(rt)){
      g1 <- gsub("\\(\\d+\\)","",rt$gene1[i])
      g2 <- gsub("\\(\\d+\\)","",rt$gene2[i])
      fuscaller <- rt$caller[i]
      num_passed <- as.integer(subset(FTSpassed,gene1==g1 & gene2==g2 & caller==fuscaller)$rec)
      if(length(num_passed) < 1) next
      total <- as.integer(subset(RSlib,gene1==g1 & gene2==g2 & caller==fuscaller)$rec)
      rt$RS[i] <- num_passed/total
    }
  return(rt)
}

# Count reads in a region using featureCounts
getBreakpointRegionCounts <- function(s){
  rt <- subset(resultTable,sample==s)
  bampath <- paste0(outputfolder,"/mapping/",s,"/",s,".bam")
  offset <- 5

  # Define regions for counting based on strand and offset
    break5regions <- cbind(paste0(rt$break5prime,"_5"),do.call(rbind.data.frame, strsplit(rt$break5prime,":")))
    colnames(break5regions) <- c("GeneID","Chr","Start","Strand")
    break5regions$Start[break5regions$Strand=="+"] <- as.numeric(break5regions$Start[break5regions$Strand=="+"]) - offset
    break5regions$Start[break5regions$Strand=="-"] <- as.numeric(break5regions$Start[break5regions$Strand=="-"]) + offset
    break5regions$End <- break5regions$Start
    
    break3regions <- cbind(paste0(rt$break3prime,"_3"),do.call(rbind.data.frame, strsplit(rt$break3prime,":")))
    colnames(break3regions) <- c("GeneID","Chr","Start","Strand")
    break3regions$Start[break3regions$Strand=="+"] <- as.numeric(break3regions$Start[break3regions$Strand=="+"]) + offset
    break3regions$Start[break3regions$Strand=="-"] <- as.numeric(break3regions$Start[break3regions$Strand=="-"]) - offset
    break3regions$End <- break3regions$Start
    
  # Perform featureCounts on region table
    regions <- unique(rbind(break5regions,break3regions))
    regions <- regions[,c(1:3,5,4)]
    temp <- with_output_sink("/dev/null", featureCounts(files = bampath, annot.ext = regions, isPairedEnd = T, ignoreDup = T, 
                                                        useMetaFeatures = F, allowMultiOverlap = T, countChimericFragments = F, 
                                                        nthreads = 2, tmpDir = resultfolder))
  
  breakpoint_counts <- list()
  idx5 <- grepl("_5",rownames(temp$counts))
  breakpoint_counts[["break5prime"]] <- temp$counts[idx5,1]
  names(breakpoint_counts[["break5prime"]]) <- gsub("_5","",rownames(temp$counts)[idx5])
  idx3 <- grepl("_3",rownames(temp$counts))
  breakpoint_counts[["break3prime"]] <- temp$counts[idx3,1]
  names(breakpoint_counts[["break3prime"]]) <- gsub("_3","",rownames(temp$counts)[idx3])
  return(breakpoint_counts)
}

calcTPM <- function(count,length,totalrpk){
  length <- length/1000
  rpk <- count/length
  scalingfactor <- totalrpk/1000000
  tpm <- rpk/scalingfactor
  return(tpm)
}

calcRPKMatrix <- function(c,gl){
  gl <- gl[match(rownames(c),gl$geneID),]
  gl$length <- gl$length/1000
  rpk <- sweep(c, MARGIN = 1, STATS = gl$length, FUN = '/')
  return(rpk)
}

calcTPMMatrix <- function(rpk){
  totalrpk <- colSums(rpk)
  totalrpk <- totalrpk/1000000
  tpm <- sweep(rpk, MARGIN = 2, STATS = totalrpk, FUN = '/')
  return(tpm)
}

zigzag_sort <- function(x, cores){
  sortvec <- rep(c(seq(1,cores),seq(cores, 1)), length = length(x))
  sortvec <- order(sortvec)
  x <- x[sortvec]
  return(x)
}
