# Helper function for calculation of the FTS
# Input:
# rt = Resulttable as produced by workflow.R
# c = raw read counts
# gl = list of gene lengths
# ism = list of median insert sizes for each sample
calcFTS <- function(rt,c,gl,ism){
  gl <- gl[match(rownames(c),gl$geneID),]
  gl$length <- gl$length/1000
  rpk <- sweep(c,MARGIN=1,gl$length,'/')
  totalrpk <- colSums(rpk)
  totalrpk <- totalrpk/1000000
  
  #normalize fusion transcripts by median insert size - TPM
  rt$fus_rpk <- 0
  rt$tpm_fusion <- 0
  for(i in 1:nrow(rt)){
    rt$fus_rpk[i] <- rt$cov[i]/(ism[ism$sample==rt$sample[i],]$insert/1000)
    rt$tpm_fusion[i] <- rt$fus_rpk[i]/totalrpk[rt$sample[i]]
  }
  rt$FTS5 <- rt$tpm_fusion/(rt$tpm_fusion+rt$tpm_gene1)
  rt$FTS3 <- rt$tpm_fusion/(rt$tpm_fusion+rt$tpm_gene2)
  rt$FTS <- round((rt$FTS5+rt$FTS3)/2,2)
  
  rt <- rt[,-which(grepl("fus_rpk",colnames(rt)))]
  return(rt)
}
