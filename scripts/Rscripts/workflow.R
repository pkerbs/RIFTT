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
                          "reshape2","ggplot2","grid","plotly",
                          "htmlwidgets","scales","future.apply",
                          "withr","Rsubread","data.table","ComplexHeatmap")
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
    debug_flag <- as.numeric(arguments[2])
    internal_BL <- as.numeric(arguments[3])
    threads <- as.numeric(arguments[4])
    
    # TODO: Implement parameter in FP_filter.sh
        cohort_flag <- T   
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
                 "hgncSymbols.R",
                 "handleCallerResults.R",
                 "importCounts.R",
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
                              cov=integer(), stringsAsFactors = F
    )
    filtered_AR <- filterResults(AR_results,"AR")
    resultTable <- rbind(resultTable,filtered_AR, stringsAsFactors = F)
    filtered_FC <- filterResults(FC_results,"FC")
    resultTable <- rbind(resultTable,filtered_FC, stringsAsFactors = F)
    
    # Define genelabels
    # In case of intergenic breakpoint:
    # Remove "(distance from flanking gene)", Convert "IGHx" names to "IGH"
      genes1 <- gsub("\\(\\d+\\)","",resultTable$gene1)
      genes1 <- gsub("^IGH.*","IGH",genes1)
      genes2 <- gsub("\\(\\d+\\)","",resultTable$gene2)
      genes2 <- gsub("^IGH.*","IGH",genes2)
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
  # Determine Promiscuity Score for every cohort individually. Else -> don't distinguish between cohorts
  if(cohort_flag){
    temp <- data.frame(cohort=character(),
                       sample=character(),
                       caller=character(),
                       gene1=character(),
                       gene2=character(),
                       break5prime=character(),
                       break3prime=character(),
                       cov=integer(),
                       label=character(),
                       known=character(),
                       reciprocal=logical(),
                       karyo=character(),
                       mol=character(),
                       PS=numeric(), stringsAsFactors = F)
    for(ch in unique(resultTable$cohort)){
        temp <- rbind(temp,addPSToTable(resultTable[resultTable$cohort==ch,]))
    }
    resultTable <- temp
  } else {
    resultTable <- addPSToTable(resultTable)
  }
  cat(" ",green("\u2713"),"\n",sep="")
#------------------------------------------------------------------------------------------------------------

# CALCULATE TPM and FTS
  cat("  ",yellow("->")," Import read counts...",sep="")
    rawcounts <- getRawCounts(paste0(outputfolder,"/featurecounts/reformatted/"))
  cat(" ",green("\u2713"),"\n",sep="")
  
  # Read in median insert sizes
  cat("  ",yellow("->")," Import estimated insert sizes...",sep="")
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
  cat(" ",green("\u2713"),"\n",sep="")
  #----------------------------

  cat("  ",yellow("->")," Calculate Fusion Transcript Scores...",sep="")
    # Run samples in parallel using future package
      run_order <- resultTable %>%
        group_by(sample) %>%
        summarise(num_pos=length(unique(c(break5prime,break3prime)))) %>%
        arrange(desc(num_pos))
      run_order <- zigzag_sort(run_order$sample, threads/2)
      
      future::plan("future::multisession", workers = threads/2)
      breakpoint_counts <- future.apply::future_lapply(run_order,
                                                       getBreakpointRegionCounts, 
                                                       future.seed= NULL)
      future::plan("sequential")
      names(breakpoint_counts) <- run_order

    # Fill resultTable with TPM values of breakpoint counts
      rpkMatrix <- calcRPKMatrix(rawcounts,genelengths)
      temp <- apply(resultTable,1,function(x){
        s <- x['sample']
        break5prime <- x['break5prime']
        break3prime <- x['break3prime']
        fusion_coverage <- as.numeric(x['cov'])
        insertsize <- inserts$insert[inserts$sample==s]
        totalrpk <- sum(rpkMatrix[,colnames(rpkMatrix)==s])
        tpm5prime <- calcTPM(breakpoint_counts[[s]]$break5prime[break5prime],insertsize,totalrpk)
        tpm3prime <- calcTPM(breakpoint_counts[[s]]$break3prime[break3prime],insertsize,totalrpk)
        tpmfusion <- calcTPM(fusion_coverage,insertsize,totalrpk)
        c(tpm5prime, tpm3prime, tpmfusion)
      })
      resultTable$tpm5prime_new <- temp[1,]
      resultTable$tpm3prime_new <- temp[2,]
      resultTable$tpmfusion_new <- temp[3,]
      resultTable$FTS5 <- resultTable$tpmfusion_new/(resultTable$tpmfusion_new+resultTable$tpm5prime_new)
      resultTable$FTS3 <- resultTable$tpmfusion_new/(resultTable$tpmfusion_new+resultTable$tpm3prime_new)
      resultTable$FTS <- (resultTable$FTS5+resultTable$FTS3)/2
   cat(" ",green("\u2713"),"\n",sep="")
#------------------------------------------------------------------------------------------------------------

# FILTERS
  # Fusion events that passed built-in filters of the callers
  # -> Evidence level 1
    resultTable$ev_level <- 1
  
  # BLACKLIST
  # Fusion events that were not detected in healthy samples
    BLpassIdx <- !resultTable$label %in% blacklist
    resultTable$ev_level[BLpassIdx] <- resultTable$ev_level[BLpassIdx] + 1
  # //DEPRECATED resultTable$ev_level[!resultTable$label %in% blacklist] <- 2

  # PS
  # Fusion events with lower PS than highest PS of known fusions among all samples
    temp <- resultTable$PS[resultTable$known=="known"]
    if(isEmpty(temp) || sum(temp) < 2){
      # If there are too few samples harboring known fusions or too few samples
      # analyzed by this filter pipeline, no reasonable PS will be determined.
      # Therefore, maxPS is adjusted to a conservative value of 2.
      maxPS <- 2
    }else{
        maxPS <- max(temp)
    }
    resultTable$ev_level[resultTable$PS <= maxPS] <- resultTable$ev_level[resultTable$PS <= maxPS] + 1
    # //DEPRECATED resultTable$ev_level[resultTable$ev_level==2
    #                                   & resultTable$PS <= maxPS] <- 3
  
  # FTS
  # Known fusions are generally regarded as highly relevant
  # and are not affected by the FTS. Revision may be required.
    # resultTable$ev_level[resultTable$known=="known"] <- resultTable$ev_level[resultTable$known=="known"]+1
    # //DEPRECATED resultTable$ev_level[resultTable$known=="known"] <- 4
  
  # Unknown fusion events are filtered by
  # 1. FTS of the 5' and 3' partner gene should be at least 0.025, which translates
  # to a relative abundance of 5% of the fusion harboring clone in a bulk RNA-seq sample
  # under the assumption of a monoallelic occurrence of the fusion.
  # 2. FTS of the 5' and 3' partner gene should be less than 1, since these events
  # suggest that there is no transcription of the single partner genes of a fusion
  # although there are fusion transcript supporting reads, which is highly unlikely.
  # 3. Mean FTS >= 0.1
    FTSpassIdx <- which(resultTable$FTS5 >= 0.025 & resultTable$FTS3 >= 0.025
                        & resultTable$FTS5 < 1 & resultTable$FTS3 < 1
                        & resultTable$FTS >= 0.05)
    resultTable$ev_level[FTSpassIdx] <- resultTable$ev_level[FTSpassIdx] + 1
    # //DEPRECATED:
    # resultTable$ev_level[resultTable$ev_level == 3
    # & resultTable$known=="unknown"
    # & resultTable$FTS5 >= 0.025 & resultTable$FTS3 >= 0.025
    # & resultTable$FTS5 < 1 & resultTable$FTS3 < 1
    # & resultTable$FTS >= 0.1
    # ] <- 4
  
  # RS
  cat("  ",yellow("->")," Calculate Robustness Scores...",sep="")
    resultTable <- addRSToTable(resultTable,FTSpassIdx)
    resultTable$ev_level[resultTable$RS >= 0.5] <- resultTable$ev_level[resultTable$RS >= 0.5] + 1
    # //DEPRECATED resultTable$ev_level[resultTable$ev_level==4 & resultTable$RS >= 0.5] <- 5
  cat(" ",green("\u2713"),"\n",sep="")
#----------------------------

# IDENTIFY OVERLAP BETWEEN ARRIBA AND FUSIONCATCHER
  cat("  ",yellow("->")," Identify overlap between Arriba and FusionCatcher calls...",sep="")
    resultTable$callerOverlap <- F
    rt_temp <- resultTable[,match(c("sample","caller","label"),colnames(resultTable))]
    resultTable$callerOverlap <- dupsBetweenGroups(rt_temp,"caller")
    rm(rt_temp)
    resultTable$ev_level[resultTable$callerOverlap == T] <- resultTable$ev_level[resultTable$callerOverlap == T] + 1
    # //DEPRECATED resultTable$ev_level[resultTable$ev_level==5 & resultTable$callerOverlap == T] <- 6
  cat(" ",green("\u2713"),"\n",sep="")
#------------------------------------------------------------------------------------------------------------

# OUTPUT RESULTS AND PLOTS
  dir.create(paste0(resultfolder,"/plots"),showWarnings = T,recursive = F)
  
  # Generate plots
  cat("  ",yellow("->")," Generate plots...",sep="")
    ps_violinplot(resultTable,BLpassIdx,paste0(resultfolder,"/plots/Violinplot_PS.png"))
    
    fts_violinplot(resultTable,paste0(resultfolder,"/plots/Violinplot_FTS.png"))
    
    saveWidget(create_TPM_FTS_plotly(resultTable),paste0(resultfolder,"/plots/TPM-FTS_3D_plot.html"))
    
    invisible(oncoprint(paste0(resultfolder,"/plots/Oncoprint_Karyo_MDx_RNAseq.png")))
    # png(paste0(resultfolder,"/plots/barplot_excluded_fusions.png"),width=7,height=9,res=600,units="cm")
    #   filterBarplot(resultTable,AR_results,FC_results)
    # invisible(dev.off()) 
    
    for(ch in unique(resultTable$cohort)){
      circosPlotsCandidates(resultTable,ch,"known")
      circosPlotsCandidates(resultTable,ch,"unknown")
    }
    
  cat(" ",green("\u2713"),"\n",sep="")
  
  cat("  ",yellow("->")," Write result table...",sep="")
    OrderIdx <- order(resultTable$known,-resultTable$ev_level)
    resultTable <- resultTable[OrderIdx,]
    # Fusion events are often reported several times (e.g. with different breakpoints)
    # Select the fusion event with the highest evidence level
    resultTable <- resultTable[!duplicated(resultTable[,c("sample","caller","label")]),]
    write.xlsx(resultTable,paste0(resultfolder,"/resultTable.xlsx"))
    # Create R workspace (only if debug flag is set)
    if(as.logical(debug_flag)){
      save.image(paste0(resultfolder,"/filterrun.RData"))
    }
  cat(" ",green("\u2713"),"\n",sep="")
#------------------------------------------------------------------------------------------------------------
