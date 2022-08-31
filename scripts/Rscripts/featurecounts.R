library(Rsubread)
library(withr)

# Initialisation
  arguments <- commandArgs(trailingOnly = T)
  sample <- arguments[1]
  annofile <- arguments[2]
  threads <- arguments[3]
  outputfolder <- "/outputfolder"
  file_fc_calls <- paste0(outputfolder,"/fusioncatcher/",sample,"/final-list_candidate-fusion-genes.txt")
  file_ar_calls <- paste0(outputfolder,"/arriba/",sample,"/fusions.tsv")
  bampath <- paste0(outputfolder,"/mapping/",sample,"/",sample,".bam")
  tempdir <- paste0(outputfolder,"/mapping/",sample,"/")
  outfile <- paste0(outputfolder,"/featurecounts/",sample,"_counts.RDS")

# Breakpoint positions are used to collect coverage near the breakpoints (required for FTS)
# Offset from the breakpoint is 10 bases
  offset <- 10
  regions <- data.frame(GeneID=character(),Chr=character(),Start=numeric(),End=numeric(),Strand=character())

  # FusionCatcher calls
    fc_calls <- read.delim(file_fc_calls, header = T, sep = "\t", check.names = F, stringsAsFactors = F)
    
    if(nrow(fc_calls) > 0){
      fc_calls$`Fusion_point_for_gene_1(5end_fusion_partner)` <- paste0("chr",fc_calls$`Fusion_point_for_gene_1(5end_fusion_partner)`)
      fc_calls$`Fusion_point_for_gene_2(3end_fusion_partner)` <- paste0("chr",fc_calls$`Fusion_point_for_gene_2(3end_fusion_partner)`)
      
      temp5 <- lapply(fc_calls$`Fusion_point_for_gene_1(5end_fusion_partner)`,function(x){
        strsplit(x,":")
      })
      temp5 <- as.data.frame(do.call(rbind,do.call(rbind, temp5)))
      colnames(temp5) <- c("Chr","Start","Strand")
      temp5$Start <- as.numeric(temp5$Start)
      temp5$Start[temp5$Strand=="+"] <- temp5$Start[temp5$Strand=="+"] - offset
      temp5$Start[temp5$Strand=="-"] <- temp5$Start[temp5$Strand=="-"] + offset
      temp5$End <- temp5$Start
      temp5$GeneID <- paste0(fc_calls$`Fusion_point_for_gene_1(5end_fusion_partner)`,"_5")
      temp5 <- temp5[,c(5,1,2,4,3)]
      
      temp3 <- lapply(fc_calls$`Fusion_point_for_gene_2(3end_fusion_partner)`,function(x){
        strsplit(x,":")
      })
      temp3 <- as.data.frame(do.call(rbind,do.call(rbind, temp3)))
      colnames(temp3) <- c("Chr","Start","Strand")
      temp3$Start <- as.numeric(temp3$Start)
      temp3$Start[temp3$Strand=="+"] <- temp3$Start[temp3$Strand=="+"] + offset
      temp3$Start[temp3$Strand=="-"] <- temp3$Start[temp3$Strand=="-"] - offset
      temp3$End <- temp3$Start
      temp3$GeneID <- paste0(fc_calls$`Fusion_point_for_gene_2(3end_fusion_partner)`,"_3")
      temp3 <- temp3[,c(5,1,2,4,3)]
      
      regions <- rbind(temp5,temp3) 
    }
  
  # Arriba calls  
    ar_calls <- read.delim(file_ar_calls, header = T, sep = "\t", check.names = F, stringsAsFactors = F)
    
    if(nrow(ar_calls) > 0){
      strand5 <- unlist(lapply(ar_calls$`strand1(gene/fusion)`,function(x){
        unlist(strsplit(x,"/"))[2]
      }))
      temp5 <- lapply(ar_calls$breakpoint1,function(x){
        strsplit(x,":")
      })
      temp5 <- as.data.frame(do.call(rbind,do.call(rbind, temp5)))
      colnames(temp5) <- c("Chr","Start")
      temp5$Start <- as.numeric(temp5$Start)
      temp5$Strand <- unlist(lapply(ar_calls$`strand1(gene/fusion)`,function(x){
        unlist(strsplit(x,"/"))[2]
      }))
      temp5$Start[temp5$Strand=="+"] <- temp5$Start[temp5$Strand=="+"] - offset
      temp5$Start[temp5$Strand=="-"] <- temp5$Start[temp5$Strand=="-"] + offset
      temp5$End <- temp5$Start
      temp5$GeneID <- paste0(ar_calls$breakpoint1,":",temp5$Strand,"_5")
      temp5 <- temp5[,c(5,1,2,4,3)]
      
      strand3 <- unlist(lapply(ar_calls$`strand2(gene/fusion)`,function(x){
        unlist(strsplit(x,"/"))[2]
      }))
      temp3 <- lapply(ar_calls$breakpoint2,function(x){
        strsplit(x,":")
      })
      temp3 <- as.data.frame(do.call(rbind,do.call(rbind, temp3)))
      colnames(temp3) <- c("Chr","Start")
      temp3$Start <- as.numeric(temp3$Start)
      temp3$Strand <- unlist(lapply(ar_calls$`strand2(gene/fusion)`,function(x){
        unlist(strsplit(x,"/"))[2]
      }))
      temp3$Start[temp3$Strand=="+"] <- temp3$Start[temp3$Strand=="+"] + offset
      temp3$Start[temp3$Strand=="-"] <- temp3$Start[temp3$Strand=="-"] - offset
      temp3$End <- temp3$Start
      temp3$GeneID <- paste0(ar_calls$breakpoint2,":",temp3$Strand,"_3")
      temp3 <- temp3[,c(5,1,2,4,3)]
      
      regions <- rbind(regions,rbind(temp5,temp3))
    }
  regions <- regions[!duplicated(regions$GeneID),]
  
  # In case there were no calls either from Arriba nor FusionCatcher
  # Stop the script here
    if(nrow(regions) == 0){
      quit(save="no")
    }
  
  temp <- with_output_sink("/dev/null", featureCounts(files = bampath, annot.ext = regions, isPairedEnd = T, countReadPairs = T, 
                                                      ignoreDup = T, useMetaFeatures = F, allowMultiOverlap = T, 
                                                      countChimericFragments = F, nthreads = threads, tmpDir = tempdir))
  
  counts <- list()
  idx5 <- grepl("_5",rownames(temp$counts))
  counts[["break5prime"]] <- temp$counts[idx5,1]
  names(counts[["break5prime"]]) <- gsub("_5","",rownames(temp$counts)[idx5])
  idx3 <- grepl("_3",rownames(temp$counts))
  counts[["break3prime"]] <- temp$counts[idx3,1]
  names(counts[["break3prime"]]) <- gsub("_3","",rownames(temp$counts)[idx3])  
  
  temp <- with_output_sink("/dev/null", featureCounts(files = bampath, annot.ext = annofile, isGTFAnnotationFile = T, 
                                                      isPairedEnd = T, countReadPairs = T, ignoreDup = T, useMetaFeatures = T, 
                                                      largestOverlap = T, primaryOnly = T, countChimericFragments = T, 
                                                      GTF.featureType = c("exon","CDS","UTR"), nthreads = threads, 
                                                      tmpDir = tempdir))
  
  counts[["allcounts"]] <- temp$counts[,1]
  names(counts[["allcounts"]]) <- rownames(temp$counts)
  
  saveRDS(counts,outfile)