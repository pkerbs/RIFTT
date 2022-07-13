# Helper function for filtering by built-in filter of fusion callers
filterResults <- function(fusionCalls_table,caller){
  minreads <- 3

  #FusionCatcher
  filterterms <- paste("(1000genomes",
                       "gap","gtex","partial","bodymap2","cacg","cortex","duplicates",
                       "ensembl_fully_overlapping","ensembl_partially_overlapping",
                       "ensembl_same_strand_overlapping","hpa","non_cancer_tissues",
                       "non_tumor_cells","pair_pseudo_genes","paralogs","refseq_fully_overlapping",
                       "refseq_partially_overlapping","refseq_same_strand_overlapping",
                       "ucsc_fully_overlapping","ucsc_partially_overlapping",
                       "ucsc_same_strand_overlapping",
                       "pseudo)",sep="|")
  if(caller=="FC"){
    keep <- c()
    for(i in 1:nrow(fusionCalls_table)){
      if(!grepl(filterterms,fusionCalls_table$Fusion_description[i]) 
         & fusionCalls_table$cov[i] >= minreads
         & fusionCalls_table$cov[i] > (2*fusionCalls_table$Counts_of_common_mapping_reads[i])){
        keep <- c(keep,i)
      }
    }
    temp <- fusionCalls_table[keep,]
  }
  #Arriba
  else if(caller=="AR"){
    temp <- subset(fusionCalls_table,confidence %in% c("medium","high") & cov>=minreads)
  }
  temp$caller <- caller
  return(temp[,match(c("cohort","sample","caller","gene1","gene2","break5prime","break3prime","cov"),colnames(temp))])
}


# import results from FusionCatcher
getFCResults <- function(folder,hgnc_symbols){
  fusioncatcher_results <- matrix(, nrow = 0, ncol = 17)
  for(s in clintable$sample){
    path_to_file <- paste0(folder,"/",s,"/final-list_candidate-fusion-genes.txt")
    tt <- read.delim(path_to_file, stringsAsFactors=F,check.names = F)
    tt <- cbind(c(rep(s,nrow(tt))),tt)
    colnames(tt)[1] <- "sample"
    colnames(tt)[which(colnames(tt)=="Gene_1_symbol(5end_fusion_partner)")] <- "gene1"
    colnames(tt)[which(colnames(tt)=="Gene_2_symbol(3end_fusion_partner)")] <- "gene2"
    colnames(tt)[which(colnames(tt)=="Gene_1_id(5end_fusion_partner)")] <- "gene_id1"
    colnames(tt)[which(colnames(tt)=="Gene_2_id(3end_fusion_partner)")] <- "gene_id2"
    fusioncatcher_results <- rbind(fusioncatcher_results,tt)
  }
  
  colnames(fusioncatcher_results)[colnames(fusioncatcher_results)=="Fusion_point_for_gene_1(5end_fusion_partner)"] <- "break5prime"
  colnames(fusioncatcher_results)[colnames(fusioncatcher_results)=="Fusion_point_for_gene_2(3end_fusion_partner)"] <- "break3prime"
  fusioncatcher_results$break5prime <- paste0("chr",fusioncatcher_results$break5prime)
  fusioncatcher_results$break3prime <- paste0("chr",fusioncatcher_results$break3prime)
  
  # DEPRECATED:
  # Separate breakpoints into columns: "chr","pos"
  #   break1 <- strsplit(fusioncatcher_results[,which(colnames(fusioncatcher_results)=="Fusion_point_for_gene_1(5end_fusion_partner)")],":")
  #   break2 <- strsplit(fusioncatcher_results[,which(colnames(fusioncatcher_results)=="Fusion_point_for_gene_2(3end_fusion_partner)")],":")
  #   fusioncatcher_results$chr1 <- unlist(lapply(break1,function(x) x[1]))
  #   fusioncatcher_results$chr2 <- unlist(lapply(break2,function(x) x[1]))
  #   fusioncatcher_results$pos1 <- unlist(lapply(break1,function(x) x[2]))
  #   fusioncatcher_results$pos2 <- unlist(lapply(break2,function(x) x[2]))
  
  # Make gene names consistent between FusionCatcher and Arriba
    temp1 <- ensemblIDToSymbol(fusioncatcher_results$gene_id1,geneNamesFromAnno)
    temp2 <- ensemblIDToSymbol(fusioncatcher_results$gene_id2,geneNamesFromAnno)
    # Keep gene name assigned by FusionCatcher if matching ENSEMBL ID was not found
      x <- which(is.na(temp1) | temp1=="")
      temp1[x] <- fusioncatcher_results$gene1[x]
      x <- which(is.na(temp2) | temp1=="")
      temp2[x] <- fusioncatcher_results$gene2[x]
    fusioncatcher_results$gene1 <- temp1
    fusioncatcher_results$gene2 <- temp2
    
  fusioncatcher_results$cov <- fusioncatcher_results$Spanning_pairs+fusioncatcher_results$Spanning_unique_reads

  return(fusioncatcher_results)
}

# import results from Arriba
getARResults <- function(folder,hgnc_symbols){
  arriba_results <- matrix(, nrow = 0, ncol = 31)
  for(s in clintable$sample){
    path_to_file <- paste0(folder,"/",s,"/fusions.tsv")
    tt <- read.delim(path_to_file, stringsAsFactors=F,check.names = F)
    tt <- cbind(c(rep(s,nrow(tt))),tt)
    colnames(tt)[1] <- "sample"
    colnames(tt)[2] <- "gene1"
    arriba_results <- rbind(arriba_results,tt)
  }
  
  # Split fusion events that have an intergenic breakpoint
  # indicated by the string: "(distance from flanking gene)" behind the genename 
    g1InterBreakidx <- which(grepl(",",arriba_results$gene1))
    flanks <- strsplit(arriba_results$gene1[g1InterBreakidx],",")
    arribaExtend1 <- arriba_results[rep(g1InterBreakidx,times=lengths(flanks)),]
    arribaExtend1$gene1 <- unlist(flanks)

    g2InterBreakidx <- which(grepl(",",arriba_results$gene2))
    flanks <- strsplit(arriba_results$gene2[g2InterBreakidx],",")
    arribaExtend2 <- arriba_results[rep(g2InterBreakidx,times=lengths(flanks)),]
    arribaExtend2$gene2 <- unlist(flanks)

    arriba_results <- arriba_results[-c(g1InterBreakidx,g2InterBreakidx),]
    arriba_results <- rbind(arriba_results,arribaExtend1)
    arriba_results <- rbind(arriba_results,arribaExtend2)
  
  # Reformat breakpoints
    colnames(arriba_results)[colnames(arriba_results)=="breakpoint1"] <- "break5prime"
    colnames(arriba_results)[colnames(arriba_results)=="breakpoint2"] <- "break3prime"
    strands5prime <- strsplit(arriba_results$`strand1(gene/fusion)`,"/")
    strands3prime <- strsplit(arriba_results$`strand2(gene/fusion)`,"/")
    strands5prime <- unlist(lapply(strands5prime,function(x) x[2]))
    strands3prime <- unlist(lapply(strands3prime,function(x) x[2]))
    arriba_results$break5prime <- paste0(arriba_results$break5prime,":",strands5prime)
    arriba_results$break3prime <- paste0(arriba_results$break3prime,":",strands3prime)
    
  # DEPRECATED 
  # Separate breakpoints into columns: "chr","pos"
  #   break1 <- strsplit(arriba_results[,which(colnames(arriba_results)=="breakpoint1")],":")
  #   break2 <- strsplit(arriba_results[,which(colnames(arriba_results)=="breakpoint2")],":")
  #   arriba_results$chr1 <- gsub("chr","",unlist(lapply(break1,function(x) x[1])))
  #   arriba_results$chr2 <- gsub("chr","",unlist(lapply(break2,function(x) x[1])))
  #   arriba_results$pos1 <- unlist(lapply(break1,function(x) x[2]))
  #   arriba_results$pos2 <- unlist(lapply(break2,function(x) x[2]))

  arriba_results$cov <- arriba_results$split_reads1+arriba_results$split_reads2+arriba_results$discordant_mates
  
  return(arriba_results)
}

# IMPORT RESULTS FROM FUSION CALLERS
cat("  ",yellow("->")," Import fusion calls...",sep="")
  # Arriba
    folder <- paste0(outputfolder,"/arriba/")
    AR_results <- getARResults(folder,hgnc_symbols)

  # FusionCatcher
    folder <- paste0(outputfolder,"/fusioncatcher/")
    FC_results <- getFCResults(folder,hgnc_symbols)

  #reduce table to the isoforms of the fusions with the highest coverage  
    # AR_results <- AR_results[order(AR_results$cov,decreasing = T),]
    # AR_results <- AR_results[!duplicated(AR_results[,c("sample","gene1","gene2")]),]
    # 
    # FC_results <- FC_results[order(FC_results$cov,decreasing = T),]
    # FC_results <- FC_results[!duplicated(FC_results[,c("sample","gene1","gene2")]),]
  
  # Assign cohort info to samples using information from clinical table
    AR_results$cohort <- clintable$cohort[match(AR_results$sample,clintable$sample)]
    FC_results$cohort <- clintable$cohort[match(FC_results$sample,clintable$sample)]
cat(" ", green("\u2713"),"\n",sep="")
#------------------------------------------------------------------------------------------------------------
