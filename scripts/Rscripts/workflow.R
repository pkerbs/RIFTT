# INITIALIZATION STEP
  library("crayon")
  cat("  ",yellow("->")," Initialize filtering pipeline...",sep="")
	# Singularity mounts the HOME directory of the user by default.
	# This leads to R using libraries installed in the HOME folder of the user. 
	# This line assures that libraries from the Singularity container are used instead.
	  .libPaths("/usr/local/lib/R/site-library")
	
	#t <- paste0(format(Sys.time(), "%H"),"_",format(Sys.time(), "%M"),"_",format(Sys.time(), "%S"))
    #d <- paste0(format(Sys.Date(), "%d"),format(Sys.Date(), "%m"),format(Sys.Date(), "%y"))
  
  # Load required packages
    list.of.packages <- c("openxlsx","dplyr","stringr",
                          "ggplot2","grid","plotly",
                          "htmlwidgets","scales","withr",
                          "data.table","ComplexHeatmap")
    # bioconductor.packages <- c("biomaRt")
    # list.of.packages <- c(list.of.packages,bioconductor.packages)
    for(p in list.of.packages){
      suppressPackageStartupMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE))
    }
    
  # Find base directory
    initial.options <- commandArgs(trailingOnly = FALSE)
    file.arg.name <- "--file="
    script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
    script.basename <- dirname(script.name)
  
  # GET ARGUMENTS
    outputfolder <- "/outputfolder"
    clintable <- read.xlsx("/clintable.xlsx") %>%
      mutate(across(everything(), as.character))
    arguments <- commandArgs(trailingOnly = T)
    resultfolder <- normalizePath(arguments[1])
    internal_BL <- as.numeric(arguments[2])
  cat(" ",green("\u2713"),"\n",sep="")
  
  # Calculate gene lengths from annotation file (printGeneLengths.jar)
  # Get EnsemblID to GeneName association table from annotation
    cat("  ",yellow("->")," Extract gene names and lengths from annotation...",sep="")
      system(paste0("java -jar /scripts/printGeneLengths.jar /anno.gtf ",resultfolder,"/genelengths.txt"))
      genelengths <- read.delim(paste0(resultfolder,"/genelengths.txt"), stringsAsFactors=FALSE)
      geneNamesFromAnno <- genelengths[,c(4,5)]
      geneNamesFromAnno$geneName <- paste0(geneNamesFromAnno$geneName,gsub("ENSG[^_]*","",geneNamesFromAnno$geneID))
    cat(" ",green("\u2713"),"\n",sep="")
  
  # Import all scripts in base directory
    scripts <- c("checkForFiles.R",
                 "importFusionDBs.R",
                 "handleCallerResults.R",
                 "calcScores.R",
                 "dupsBetweenGroups.R",
                 "plots.R"
                 )
    for(s in scripts){
      script <- file.path(script.basename, s)
      source(script)
    }
  # GENERATE BLACKLIST FROM FUSIONS FOUND IN HEALTHY SAMPLES
    cat("  ",yellow("->")," Load blacklists of fusion genes...",sep="")
      blacklist <- c()
      if(as.logical(internal_BL)){
        blacklist <- read.xlsx(paste0(script.basename,"/tables/fusion_blacklist.xlsx"))
        blacklist <- blacklist[,1]
      }
      if(file.size("/blacklist_user.xlsx") > 0){
        blacklist_user <- read.xlsx("/blacklist_user.xlsx",)[,1]
        blacklist <- c(blacklist,blacklist_user)
      }
    cat(" ",green("\u2713"),"\n",sep="")
  #------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------

# FILTER RESULTS FROM FUSION CALLERS
  cat("  ",yellow("->")," Apply built-in filters...",sep="")
    resultTable <- data.frame(cohort=character(),
                              sample=character(),
                              caller=character(),
                              gene1=character(),
                              gene2=character(),
                              break5prime=character(),
                              break3prime=character(),
                              cov=integer(),
                              passCallFilt=logical(), stringsAsFactors = F
    )
    AR_results <- applyCallerFilter(AR_results,"AR")
    FC_results <- applyCallerFilter(FC_results,"FC")
    temp1 <- AR_results[,match(c("cohort","sample","caller","gene1","gene2","break5prime","break3prime","cov","passCallFilt"),colnames(AR_results))]
    temp2 <- FC_results[,match(c("cohort","sample","caller","gene1","gene2","break5prime","break3prime","cov","passCallFilt"),colnames(FC_results))]
    resultTable <- rbind(temp1,temp2)
    rm(temp1,temp2)
    
    # Define genelabels
    # In case of intergenic breakpoint:
    # Remove "(distance from flanking gene)", Convert "IGHx" names to "IGH", Convert "TR[ABCD]x" names to "TR[ABCD]"
      genes1 <- gsub("\\(\\d+\\)","",resultTable$gene1)
      genes1 <- gsub("^IGH.*","IGH",genes1)
      genes1 <- gsub("(^TR[ABCD]).*","\\1",genes1)
      genes2 <- gsub("\\(\\d+\\)","",resultTable$gene2)
      genes2 <- gsub("^IGH.*","IGH",genes2)
      genes2 <- gsub("(^TR[ABCD]).*","\\1",genes2)
      resultTable$label <- paste0(genes1,"::",genes2)
      resultTable$reciproc_label <- paste0(genes2,"::",genes1)
  cat(" ",green("\u2713"),"\n",sep="")
#------------------------------------------------------------------------------------------------------------

# INSERT INFO ABOUT KNOWN FUSION  
  cat("  ",yellow("->")," Identify known fusions and reciprocal events...",sep="")
  resultTable$known <- "unknown"
  chimerdb_known <- c(subset(chimerkb,Fusion_pair %in% chimerpub$Fusion_pair)$Fusion_pair,
                      subset(chimerpub,rec>1)$Fusion_pair)
  resultTable$known[resultTable$label %in% chimerdb_known
                    | resultTable$reciproc_label %in% chimerdb_known] <- "known"
  
  # Adjust fusion label if reciprocal version has more entries in MitelmanDB
  # Determine reciprocal events for fusions with no entry in MitelmanDB
  # If fusion and reciprocal fusion not found in MitelmanDB,
  # first compare by coverage then by lexicographical order
    resultTable$reciprocal <- F
    for(i in 1:nrow(resultTable)){
      s <- resultTable$sample[i]
      c <- resultTable$caller[i]
      f <- resultTable$label[i]
      fr <- resultTable$reciproc_label[i]
      if(f==fr) next
      if(sum(mitelman$rec[mitelman$fusion==fr]) > sum(mitelman$rec[mitelman$fusion==f])){
        if(fr %in% subset(resultTable,sample==s & caller==c)$label){
          resultTable$reciprocal[i] <- T
        }
        resultTable$label[i] <- fr
        resultTable$reciproc_label[i] <- f
        next
      }
      if(!f %in% mitelman$fusion & !fr %in% mitelman$fusion){
        if(fr %in% subset(resultTable,sample==s & caller==c)$label){
          f_cov <- resultTable$cov[i]
          fr_cov <- max(subset(resultTable,sample==s & caller==c & label==fr)$cov)
          if(fr_cov > f_cov){
            resultTable$reciprocal[i] <- T
            resultTable$label[i] <- fr
            resultTable$reciproc_label[i] <- f
          } else {
            if(fr_cov == f_cov){
              if(fr < f){
                resultTable$reciprocal[i] <- T
                resultTable$label[i] <- fr
                resultTable$reciproc_label[i] <- f
              }
            }
          }
        }
      }
    }
    resultTable$mitelman_rec <- mitelman$rec[match(resultTable$label,mitelman$fusion)]
    resultTable$mitelman_rec[is.na(resultTable$mitelman_rec)] <- 0
  cat(" ",green("\u2713"),"\n",sep="")
      
  #check if fusion was confirmed by routine
    cat("  ",yellow("->")," Compare clinical data to fusion calls...",sep="")
    resultTable$karyo <- NA
    resultTable$mol <- NA
    for(i in 1:nrow(resultTable)){
      s <- resultTable$sample[i]
      if(!s %in% clintable$sample) next
      
      f <- resultTable$label[i]
      fr <- resultTable$reciproc_label[i]
      g1 <- unlist(strsplit(f,"::"))[1]
      g2 <- unlist(strsplit(f,"::"))[2]
      otherCyto <- clintable$otherCyto[which(clintable$sample==s)]
      Karyotype <- clintable$Karyotype[which(clintable$sample==s)]
      
      # check mol
        if(grepl(f,otherCyto)
           | grepl(fr,otherCyto)){
          resultTable$mol[i] <- "Y"
        } else if(grepl(paste0(g1,"#"),otherCyto)
                   | grepl(paste0(g2,"#"),otherCyto)){
          resultTable$mol[i] <- "S"
        } else if(grepl(paste0(g1,"!",g2),otherCyto)
                  | grepl(paste0(g2,"!",g1),otherCyto)){
          resultTable$mol[i] <- "N"
        }
        
      # check Karyotype
        if(is.na(Karyotype)) next
        resultTable$karyo[i] <- "N"
        f_rec <- sum(mitelman$rec[mitelman$fusion==f])
        fr_rec <- sum(mitelman$rec[mitelman$fusion==fr])
        if(fr_rec > f_rec){
          karyogrep <- mitelman$karyogrep[mitelman$fusion==fr]
          karyogrep_sloppy <- mitelman$karyogrep_sloppy[mitelman$fusion==fr]
        } else{
          karyogrep <- mitelman$karyogrep[mitelman$fusion==f]
          karyogrep_sloppy <- mitelman$karyogrep_sloppy[mitelman$fusion==f]
        }
        
        if(length(karyogrep) > 0 && grepl(karyogrep,Karyotype,perl=T,ignore.case=T)){
          resultTable$karyo[i] <- "Y"
        }
        if(resultTable$karyo[i] != "Y"){
          if(length(karyogrep_sloppy) > 0 && grepl(karyogrep_sloppy,Karyotype,perl=T,ignore.case=T)){
            resultTable$karyo[i] <- "S"
          }
        }
    }
    cat(" ",green("\u2713"),"\n",sep="")
#------------------------------------------------------------------------------------------------------------

# CALCULATE PROMISCUITY SCORE
  cat("  ",yellow("->")," Calculate Promiscuity Scores...",sep="")
    promiscuity <- calcPromiscuity(resultTable)
    resultTable$PS <- unlist(apply(resultTable,1,function(x){
      gene1_temp <- (promiscuity[[x["gene1"]]]["AR_promiscuity"] + promiscuity[[x["gene1"]]]["FC_promiscuity"])/2
      gene2_temp <- (promiscuity[[x["gene2"]]]["AR_promiscuity"] + promiscuity[[x["gene2"]]]["FC_promiscuity"])/2
      PS <- (gene1_temp + gene2_temp)/2
      return(PS)
    }))
  cat(" ",green("\u2713"),"\n",sep="")
#------------------------------------------------------------------------------------------------------------

# CALCULATE FUSION TRANSCRIPT SCORE
  cat("  ",yellow("->")," Calculate Fusion Transcript Scores...",sep="")
    # Read in median insert sizes
      inserts <- data.frame(sample=character(),insert=numeric(),stringsAsFactors = F)
      for(s in clintable$sample){
        conn <- file(paste0(outputfolder,"/insertsizes/",s,"_is_metrics.txt"), "r")
        i <- 1
        while(length(line <- readLines(conn, 1)) > 0) {
          if(i==8) break
          i <- i+1
        }
        close(conn)
        insert <- unlist(strsplit(line,"\t"))[1]
        inserts[nrow(inserts)+1,] <- data.frame(sample=s,insert=insert,stringsAsFactors = F)
        inserts$insert <- as.numeric(inserts$insert)
      }
  
    # Import counts
      readcounts <- list()
      for(s in clintable$sample){
        readcounts[[s]] <- readRDS(paste0(outputfolder,"/featurecounts/",s,"_counts.RDS"))
      }

    # Fill resultTable with TPM values of breakpoint counts
      rpkMatrix <- calcRPKMatrix(readcounts,genelengths)
      temp <- apply(resultTable,1,function(x){
        s <- x['sample']
        break5prime <- x['break5prime']
        break3prime <- x['break3prime']
        fusion_coverage <- as.numeric(x['cov'])
        insertsize <- inserts$insert[inserts$sample==s]
        totalrpk <- sum(rpkMatrix[,colnames(rpkMatrix)==s])
        tpm5prime <- calcTPM(readcounts[[s]]$break5prime[break5prime],insertsize,totalrpk)
        tpm3prime <- calcTPM(readcounts[[s]]$break3prime[break3prime],insertsize,totalrpk)
        tpmfusion <- calcTPM(fusion_coverage,insertsize,totalrpk)
        c(tpm5prime, tpm3prime, tpmfusion)
      })
      resultTable$tpm5prime <- temp[1,]
      resultTable$tpm3prime <- temp[2,]
      resultTable$tpmfusion <- temp[3,]
      resultTable$FTS5 <- resultTable$tpmfusion/(resultTable$tpmfusion+resultTable$tpm5prime)
      resultTable$FTS3 <- resultTable$tpmfusion/(resultTable$tpmfusion+resultTable$tpm3prime)
      resultTable$FTS <- (resultTable$FTS5+resultTable$FTS3)/2
   cat(" ",green("\u2713"),"\n",sep="")
#------------------------------------------------------------------------------------------------------------

# FILTERS
   # Fusion events that passed built-in filters of the callers
   # -> Evidence level 1
   resultTable$ev_level <- 0
   resultTable$ev_level[resultTable$passCallFilt] <- resultTable$ev_level[resultTable$passCallFilt] + 1
   
   # BLACKLIST
   # Fusion events that were not detected in healthy samples
   BLpassIdx <- !paste0(resultTable$gene1,"::",resultTable$gene2) %in% blacklist
   resultTable$ev_level[BLpassIdx] <- resultTable$ev_level[BLpassIdx] + 1
   
   # PS
   # Fusion events with lower PS than highest PS of known fusions among all samples
   # If there are too few samples harboring known fusions or too few samples
   # analyzed by this filter pipeline, no reasonable PS will be determined.
   # Therefore, maxPS is adjusted to a conservative value of 2.
   temp <- resultTable$PS[resultTable$known=="known"]
   if(length(temp)==0 || sum(temp) < 2){
     maxPS <- 2
     }else{
       maxPS <- max(temp)
       }
   PSpassIdx <- resultTable$PS <= maxPS
   resultTable$ev_level[PSpassIdx] <- resultTable$ev_level[PSpassIdx] + 1

  # FTS
  # 1. FTS of the 5' and 3' partner gene should be at least 0.025, which translates
  # to a relative abundance of 5% of the fusion harboring clone in a bulk RNA-seq sample
  # under the assumption of a monoallelic occurrence of the fusion.
  # 2. FTS of the 5' and 3' partner gene should be less than 1, since these events
  # suggest that there is no transcription of the single partner genes of a fusion
  # although there are fusion transcript supporting reads, which is highly unlikely.
    FTSpassIdx <- which(resultTable$FTS5 >= 0.025 & resultTable$FTS3 >= 0.025
                        & resultTable$FTS5 < 1 & resultTable$FTS3 < 1
                        & resultTable$FTS >= 0.05)
    resultTable$ev_level[FTSpassIdx] <- resultTable$ev_level[FTSpassIdx] + 1

  # RS
  cat("  ",yellow("->")," Calculate Robustness Scores...",sep="")
    resultTable <- addRSToTable(resultTable,FTSpassIdx)
    RSpassIdx <- resultTable$RS >= 0.5
    resultTable$ev_level[RSpassIdx] <- resultTable$ev_level[RSpassIdx] + 1
  cat(" ",green("\u2713"),"\n",sep="")
#----------------------------

# IDENTIFY OVERLAP BETWEEN ARRIBA AND FUSIONCATCHER
  cat("  ",yellow("->")," Identify overlap between Arriba and FusionCatcher calls...",sep="")
    resultTable$callerOverlap <- F
    rt_temp <- resultTable[,match(c("sample","caller","label"),colnames(resultTable))]
    resultTable$callerOverlap <- dupsBetweenGroups(rt_temp,"caller")
    rm(rt_temp)
    resultTable$ev_level[resultTable$callerOverlap == T] <- resultTable$ev_level[resultTable$callerOverlap == T] + 1
  cat(" ",green("\u2713"),"\n",sep="")
#------------------------------------------------------------------------------------------------------------

# OUTPUT RESULTS AND PLOTS
  cat("  ",yellow("->")," Write result table...",sep="")
  resultTable$passBL <- F
  resultTable$passBL[BLpassIdx] <- T
  resultTable$passPS <- F
  resultTable$passPS[PSpassIdx] <- T
  resultTable$passFTS <- F
  resultTable$passFTS[FTSpassIdx] <- T
  resultTable$passRS <- F
  resultTable$passRS[RSpassIdx] <- T
  colOrder <- c("cohort","sample","caller","gene1","gene2","label","reciprocal",
                "break5prime","break3prime","cov","known","mitelman_rec","PS",
                "tpm5prime","tpm3prime","tpmfusion","FTS5","FTS3","FTS","RS","ev_level",
                "passCallFilt","passBL","passPS","passFTS","passRS","callerOverlap","karyo","mol")
  colOrderIdx <- match(colOrder,colnames(resultTable))
  resultTable <- resultTable[,colOrderIdx]
  rowOrderIdx <- order(resultTable$known,-resultTable$ev_level)
  resultTable <- resultTable[rowOrderIdx,]
  
  # Fusion events are often reported several times (e.g. with different breakpoints)
  # Select the fusion event with the highest evidence level
  resultTable <- resultTable[!duplicated(resultTable[,c("sample","caller","label")]),]
  write.xlsx(resultTable,paste0(resultfolder,"/resultTable.xlsx"))
  # Create R workspace
  save.image(paste0(resultfolder,"/filterrun.RData"))
  
  cat(" ",green("\u2713"),"\n",sep="")
  
  # Generate plots
  dir.create(paste0(resultfolder,"/plots"),showWarnings = T,recursive = F)
  cat("  ",yellow("->")," Generate plots...",sep="")
    ps_violinplot(resultTable,paste0(resultfolder,"/plots/Violinplot_PS.png"))
    
    fts_violinplot(resultTable,paste0(resultfolder,"/plots/Violinplot_FTS.png"))
    
    saveWidget(create_TPM_FTS_plotly(resultTable),paste0(resultfolder,"/plots/TPM-FTS_3D_plot.html"))
    
    invisible(oncoprint(paste0(resultfolder,"/plots/Oncoprint_Karyo_MDx_RNAseq.png")))

    for(ch in unique(resultTable$cohort)){
      rt <- subset(resultTable,known=="known"
                   & ev_level %in% c(3:6)
                   & cohort==ch
                   & reciprocal==F)
      if(nrow(rt)==0){
        cat("\n    ",yellow("->")," There are no known fusions in the '",ch,"' cohort. Circos plot skipped.",sep="")
      } else{circosPlotsCandidates(rt)}
      
      rt <- subset(resultTable,known=="unknown"
                   & ev_level %in% c(6)
                   & cohort==ch
                   & reciprocal==F)
      if(nrow(rt)==0){
        cat("\n    ",yellow("->")," There are no robust fusion candidates in the '",ch,"' cohort. Circos plot skipped.",sep="")
      } else{circosPlotsCandidates(rt)}
    }
    
  cat(" ",green("\u2713"),"\n",sep="")
#------------------------------------------------------------------------------------------------------------
