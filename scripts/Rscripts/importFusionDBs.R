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
  # mbca$karyo[mbca$fusion=="KAT6A::CREBBP"] <- paste0(mbca$karyo[mbca$fusion=="KAT6A::CREBBP"],
  #                                                   "|t(8;16)(p11.2;p13.3)")
  # mbca$karyo[mbca$fusion=="KMT2A::MLLT1"] <- gsub("p13","p13.3",mbca$karyo[mbca$fusion=="KMT2A::MLLT1"])
  # mbca$karyo[mbca$fusion=="KMT2A::ELL"] <- gsub("p13","p13.1",mbca$karyo[mbca$fusion=="KMT2A::ELL"])
  # mbca$karyo[mbca$fusion=="RPN1::MECOM"] <- paste0(mbca$karyo[mbca$fusion=="RPN1::MECOM"],
  #                                                 "|t(3;3)(q21.3;q26.2)|inv(3)(q21.3q26.2)")
  # mbca$karyo[mbca$fusion=="CBFB::MYH11"] <- paste0(mbca$karyo[mbca$fusion=="CBFB::MYH11"],
  #                                                 "|t(16;16)(p13.1;q22)|inv(16)(p13.1q22)")
  # mbca$karyo[mbca$fusion=="BCR::FGFR1"] <- paste0(mbca$karyo[mbca$fusion=="BCR::FGFR1"],
  #                                                     "|t(8;22)(p11.2;q11.2)")
  mbca$karyo[mbca$fusion=="FGFR1OP2::FGFR1"] <- paste0(mbca$karyo[mbca$fusion=="FGFR1OP2::FGFR1"],
                                                       "|t(8;12)(p12;p11)")
  mbca$karyo[mbca$fusion=="PML::RARA"] <- paste0(mbca$karyo[mbca$fusion=="PML::RARA"],
                                                  "|t(15;17)(q22;q21)")
  
  
  exceptions <- c("KMT2A::MLLT1","KMT2A::ELL")
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
  chimerpub <- subset(chimerpub,grepl("PCR|Sanger",Validation,ignore.case = T))
  chimerpub <- chimerpub %>%
    group_by(Fusion_pair) %>%
    mutate(rec=length(unique(id)))
  chimerkb <- read.xlsx(paste0(script.basename,"/tables/ChimerKB4.xlsx"))
  chimerkb <- subset(chimerkb, grepl("Plus|Curation",Source,ignore.case = T))
  chimerkb <- chimerkb %>%
    group_by(Fusion_pair) %>%
    mutate(rec=length(unique(id)))

  # FUSION LABEL FIX
    chimerkb$Fusion_pair <- paste0(chimerkb$H_gene,"::",chimerkb$T_gene)
    chimerpub$Fusion_pair <- paste0(chimerpub$H_gene,"::",chimerpub$T_gene)
  # chimerseq <- read.xlsx(paste0(script.basename,"/tables/ChimerSeq4.xlsx"))
  # chimerseq <- subset(chimerseq,Highly_Reliable_Seq=="Seq+")
cat(" ", green("\u2713"),"\n",sep="")