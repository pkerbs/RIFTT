calcPromiscuity <- function(rt){
  genes <- unique(c(rt$gene1,rt$gene2))
  genes <- genes[!duplicated(genes)]
  
  # For every gene count number of different partners
  # Count only once if different partners reported but with the same breakpoint
  promiscuity <- lapply(genes,function(x){
    temp1 <- subset(rt,caller=="AR" & gene1==x)[,grepl("gene2|break3prime",colnames(rt))]
    temp2 <- subset(rt,caller=="AR" & gene2==x)[,grepl("gene1|break5prime",colnames(rt))]
    colnames(temp1) <- c("partner","breakpoint")
    colnames(temp2) <- c("partner","breakpoint")
    temp_merged <- rbind(temp1,temp2)
    temp_merged <- temp_merged[!duplicated(temp_merged$partner),]
    temp_merged <- temp_merged[!duplicated(temp_merged$breakpoint),]
    ar_partner_num <- nrow(temp_merged)
    
    temp1 <- subset(rt,caller=="FC" & gene1==x)[,grepl("gene2|break3prime",colnames(rt))]
    temp2 <- subset(rt,caller=="FC" & gene2==x)[,grepl("gene1|break5prime",colnames(rt))]
    colnames(temp1) <- c("partner","breakpoint")
    colnames(temp2) <- c("partner","breakpoint")
    temp_merged <- rbind(temp1,temp2)
    temp_merged <- temp_merged[!duplicated(temp_merged$partner),]
    temp_merged <- temp_merged[!duplicated(temp_merged$breakpoint),]
    fc_partner_num <- nrow(temp_merged)
    
    result <- c(ar_partner_num,fc_partner_num)
    names(result) <- c("AR_promiscuity","FC_promiscuity")
    return(result)
  })
  names(promiscuity) <- genes
  return(promiscuity)
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

calcTPM <- function(count,length,totalrpk){
  length <- length/1000
  rpk <- count/length
  scalingfactor <- totalrpk/1000000
  tpm <- rpk/scalingfactor
  return(tpm)
}

calcRPKMatrix <- function(readcounts,gl){
  c <- do.call(cbind,lapply(readcounts,function(x){
    x$allcounts
  }))
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