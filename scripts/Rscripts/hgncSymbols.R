# Retrieve up-to-date EnsemblID to HGNC symbols conversion table
# getHGNCsymbols <- function(server){
#   require(biomaRt)
#   httr::set_config(httr::config(ssl_verifypeer = FALSE))
#   ensembl = useMart(host = server,
#                     biomart = "ENSEMBL_MART_ENSEMBL",
#                     dataset = "hsapiens_gene_ensembl"
#                     )
#   ens_symbols = getBM(attributes=c("ensembl_gene_id","hgnc_symbol"),mart = ensembl)
#   ens_symbols <- ens_symbols[!duplicated(ens_symbols$ensembl_gene_id),]
#   return(ens_symbols)
# }

# ENSEMBL ID to HGNC Symbol conversion for GENCODE annotation
ensemblIDToSymbol <- function(l,ens_symbols){
  tmpEnsIds <- gsub("\\.\\d+","",l)
  ens_symbols[,1] <- gsub("\\.\\d+","",ens_symbols[,1])
  gene_match <- match(tmpEnsIds,ens_symbols[,1])
  hgncIds <- ens_symbols[gene_match,2]
  return(hgncIds)
}

# cat("  ",yellow("->")," Download ENSEMBL/HGNC conversion table...",sep="")
#   hgnc_symbols <- getHGNCsymbols("https://www.ensembl.org")
# cat(" ", green("\u2713"),"\n",sep="")
