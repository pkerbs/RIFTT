reformatMitelman <- function(folder){
  require(gsubfn)
  mbca <- read.delim(paste0(folder,"/MBCA.TXT.DATA"))
  colnames(mbca)[grepl("GeneShort",colnames(mbca))] <- "fusion"
  mbca <- mbca[grepl("::",mbca$fusion),]
  mbca <- mbca[!grepl(",",mbca$fusion)
               & !grepl(",",mbca$KaryShort)
               & grepl("\\(",mbca$KaryShort),]
  mbca <- mbca %>%
    group_by(fusion) %>%
    summarise(karyo=paste0(unique(KaryShort),collapse = "|"),
              rec=length(unique(RefNo)))
  
  # Add old ISCNs for known fusions into Mitelman
    mbca$karyo[mbca$fusion=="FGFR1OP2::FGFR1"] <- paste0(mbca$karyo[mbca$fusion=="FGFR1OP2::FGFR1"],
                                                         "|t(8;12)(p12;p11)")
    mbca$karyo[mbca$fusion=="PML::RARA"] <- paste0(mbca$karyo[mbca$fusion=="PML::RARA"],
                                                    "|t(15;17)(q22;q21)")
    mbca$karyo[mbca$fusion=="PICALM::MLLT10"] <- paste0(mbca$karyo[mbca$fusion=="PICALM::MLLT10"],
                                                   "|t(10;11)(p13;q21)")
    mbca$karyo[mbca$fusion=="KMT2A::MLLT3"] <- paste0(mbca$karyo[mbca$fusion=="KMT2A::MLLT3"],
                                                        "|t(9;11)(p22;q23)")
  
  
  # exceptions <- c("KMT2A::MLLT1","KMT2A::ELL")
  mbca$karyogrep <- gsub("\\(","\\\\(",mbca$karyo)
  mbca$karyogrep <- gsub("\\)","\\\\)",mbca$karyogrep)
  mbca$karyogrep <- gsub("\\+","\\\\+",mbca$karyogrep)
  mbca$karyogrep <- gsub("([pq]\\d+)","\\1(\\\\.\\\\d+)\\?",mbca$karyogrep)
  # mbca$karyogrep[!mbca$fusion %in% exceptions] <- gsubfn("\\.\\d+",~paste0("(",x,")?"),mbca$karyogrep[!mbca$fusion %in% exceptions])
  mbca$karyo_sloppy <- gsub("\\(([pq]\\d+(\\.\\d+)?;?)+\\)","",mbca$karyo)
  mbca$karyogrep_sloppy <- gsub("\\(","\\\\(",mbca$karyo_sloppy)
  mbca$karyogrep_sloppy <- gsub("\\)","\\\\)",mbca$karyogrep_sloppy)
  mbca$karyogrep_sloppy <- gsub("\\+","\\\\+",mbca$karyogrep_sloppy)
  #--------------------------  
  return(mbca)
}

cat("  ",yellow("->")," Import ChimerDB and MitelmanDB...",sep="")
  mitCache <-  file.path(script.basename,"tables/mitelman.txt")
  mitelman <- read.delim(mitCache,check.names = F,
                         stringsAsFactors = F)
  chimerpub <- read.xlsx(paste0(script.basename,"/tables/ChimerPub4.xlsx"))
  chimerpub$Fusion_pair <- gsub("-","::",chimerpub$Fusion_pair)
  chimerpub <- chimerpub %>%
    group_by(Fusion_pair) %>%
    summarize(rec=length(unique(PMID)), val=paste0(Validation, collapse = ","))
  chimerpub$val <- unlist(lapply(chimerpub$val, function(x){
    temp <- unlist(strsplit(x,","))
    temp <- gsub("\\s", "", temp)  
    temp <- temp[temp!="NA"]
    temp <- unique(temp)
    paste0(temp,collapse = ",")
  }))
  
  chimerkb <- read.xlsx(paste0(script.basename,"/tables/ChimerKB4.xlsx"))
  chimerkb$Fusion_pair <- gsub("-","::",chimerkb$Fusion_pair)
  chimerkb <- chimerkb %>%
    group_by(Fusion_pair) %>%
    summarize(rec=paste0(PMID,collapse = ","), val=paste0(Validation, collapse = ","))
  chimerkb$rec <- unlist(lapply(chimerkb$rec, function(x){
    temp <- unlist(strsplit(x,","))
    temp <- gsub("\\s", "", temp)  
    temp <- temp[temp!="NA"]
    temp <- unique(temp)
    length(temp)
  }))
  chimerkb$val <- unlist(lapply(chimerkb$val, function(x){
    temp <- unlist(strsplit(x,","))
    temp <- gsub("\\s", "", temp)  
    temp <- temp[temp!="NA"]
    temp <- unique(temp)
    paste0(temp,collapse = ",")
  }))
cat(" ", green("\u2713"),"\n",sep="")