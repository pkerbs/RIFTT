# Promiscuity Score violin plot
ps_violinplot <- function(rt,outfile){
  require(ggplot2)
  marg <- margin(2,2,2,2,"mm")
  y_breaks <- c(0,1,2,4,6,8,10,seq(20,100,20),seq(200,1000,100))
  
  tfp <- rt[rt$passBL,]
  tfp <- tfp[!duplicated(tfp[,c("sample","label")]),]
  num_cohorts <- length(unique(tfp$cohort))
  
  p <- ggplot(tfp,aes(x=known,y=PS,fill=known)) +
    geom_violin(lwd=0.2,alpha=0.5) +
    geom_boxplot(lwd=0.2,width=0.1,outlier.shape = NA,
                 alpha=0.5, show.legend = F) +
    theme_minimal() +
    facet_wrap(~cohort,nrow = 1,strip.position = "bottom") +
    theme(legend.position = c(0.2/(num_cohorts-0.999),0.97),
          legend.key.size = unit(0.5,"line"),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 6, face="bold"),
          legend.spacing.y = unit(0.1,"mm"),
          axis.ticks.y = element_line(),
          axis.line.y = element_line(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank()
    ) +
    scale_fill_manual(values=c("#67a9cf","#ef8a62")) +
    scale_y_continuous(expand = c(0,0), trans = "log10",
                       breaks = y_breaks) +
    labs(y = "Promiscuity Score", fill="Fusion events")
  p <- p + theme(plot.margin = marg)
  
  ggsave(p, filename = outfile, 
         height=4, width=1+(1.6*num_cohorts), units = "in", dpi = 300)
}
#-----------------------------------

# Fusion Transcript Score violin plot
fts_violinplot <- function(rt,outfile){
  require(ggplot2)
  marg <- margin(2,2,2,2,"mm")
  
  tfp <- rt[which(rt$passFTS & rt$passBL),]

  p <- ggplot(tfp,aes(known,FTS,fill=known)) +
    geom_violin(lwd=0.2,alpha=0.5) +
    geom_boxplot(lwd=0.2,width=0.1,outlier.shape = NA,
                 alpha=0.5, show.legend = F) +
    theme_minimal() +
    theme(legend.position = c(0.2,0.97),
          legend.key.size = unit(0.5,"line"),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 6, face="bold"),
          legend.spacing.y = unit(0.1,"mm"),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.line.y = element_line(),
          axis.ticks.y = element_line()
    ) +
    scale_fill_manual(values = c("#67a9cf", "#ef8a62"))+
    scale_y_continuous(breaks = pretty(0:1,10),expand = c(0,0)) +
    labs(y = "Fusion Transcript Score", fill="Fusion events") 
  p <- p + theme(plot.margin = marg)
  
  ggsave(p, filename = outfile, 
         height=4, width=3, units = "in", dpi = 300)
}
#-----------------------------------

# Barplot of filtered events
filterBarplot <- function(rt,AR_results,FC_results){
  require(ggplot2)
  
  tfp <- data.frame(fgroup=character(),caller=factor(levels = c("AR","FC")),cohort=character(),
                     count=numeric(),frac=numeric(),stringsAsFactors = F)
  events <- list() 
  events[["AR"]] <- data.frame(total=0,builtin=0,healthy_filter=0,ps_filter=0,
                               fts_filter=0,rs_filter=0,candidates=0)
  events[["FC"]]<- data.frame(total=0,builtin=0,healthy_filter=0,ps_filter=0,
                              fts_filter=0,rs_filter=0,candidates=0)
  cohorts <- unique(rt$cohort)
  
  for(ch in cohorts){
    for(caller in levels(tfp$caller)){
      if(caller=="AR") total <- nrow(AR_results[AR_results$cohort==ch,])
      else total <- nrow(FC_results[FC_results$cohort==ch,])
      
      builtin_filter <- total-nrow(rt[rt$cohort==ch & rt$caller==caller,])
      healthy_filter <- nrow(rt[rt$cohort==ch & rt$caller==caller & rt$ev_level == 1,])
      ps_filter <- nrow(rt[rt$cohort==ch & rt$caller==caller & rt$ev_level == 2,])
      fts_filter <- nrow(rt[rt$cohort==ch & rt$caller==caller & rt$ev_level == 3,])
      rs_filter <- nrow(rt[rt$cohort==ch & rt$caller==caller & rt$ev_level == 4,])
      candidates <- nrow(rt[rt$cohort==ch & rt$caller==caller & rt$ev_level >= 5,])
      
      events[[caller]] <- events[[caller]] + c(total,builtin_filter,healthy_filter,ps_filter,
                                               fts_filter,rs_filter,candidates)
      tfp[nrow(tfp)+1,] <- c("Built-in Filter",caller,ch,builtin_filter,builtin_filter/total)
      tfp[nrow(tfp)+1,] <- c("Blacklist",caller,ch,healthy_filter,healthy_filter/total)
      tfp[nrow(tfp)+1,] <- c("PS Filter",caller,ch,ps_filter,ps_filter/total)
      tfp[nrow(tfp)+1,] <- c("FTS Filter",caller,ch,fts_filter,fts_filter/total)
      tfp[nrow(tfp)+1,] <- c("RS Filter",caller,ch,rs_filter,rs_filter/total)
      tfp[nrow(tfp)+1,] <- c("Candidates",caller,ch,candidates,candidates/total)
    }
  }
  tfp$count <- as.numeric(tfp$count)
  tfp$frac <- as.numeric(tfp$frac)
  tfp$fgroup <- factor(tfp$fgroup,levels=c("Built-in Filter",
                                           "Blacklist",
                                           "PS Filter",
                                           "FTS Filter",
                                           "RS Filter",
                                           "Candidates"))
  tfp$cohort <- factor(tfp$cohort,levels=cohorts)
  tfp$label <- round(100*tfp$frac,digits = 0)
  
  # calculate midpoints of bars (simplified using comment by @DWin)
  tfp <- plyr::ddply(tfp, c("cohort","caller"),
                      transform, pos = 100 - (100*((cumsum(frac) - (0.5 * frac))))
  )
  # earthcol <- c("#A16928","#bd925a","#d6bd8d","#edeac2","#b5c8b8","#79a7ac","#2887a1")
  geysircol <- c("#008080","#70a494","#b4c8a8","#f6edbd","#edbb8a","#de8a5a","#ca562c")
  tfp$caller<- gsub("AR","Arriba",tfp$caller)
  tfp$caller<- gsub("FC","FusionCatcher",tfp$caller)
  
  p<-ggplot(tfp,aes(x=caller,y=100*frac,fill=fgroup)) +
    geom_bar(stat="identity",position="stack",width=0.5,size=0.1) +
    geom_text(aes(y=pos,label=ifelse(frac<0.02,"",label)),
              size=2.8,fontface="bold") +
    facet_wrap(~cohort, strip.position = "top", nrow = 1) +
    theme_minimal() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(size=8),
          panel.grid = element_blank(),
          axis.ticks.y = element_line(),
          strip.text = element_text(size = 10,face="bold"),
          strip.background = element_rect(fill="white",linetype = "blank") ,
          legend.text = element_text(size=8),
          legend.title = element_text(size=8,face="bold"),
          legend.key.size = unit(0.3, "cm"),
          legend.key.width = unit(0.3,"cm"),
          legend.margin = margin(0,0,0,0),
          legend.box.margin = margin(0,0,-10,-10)) +
    scale_x_discrete(expand = c(0.30,0.30)) +
    scale_y_continuous(expand = c(0.01,0.01), breaks = seq(0,100,25)) +
    scale_fill_discrete(type = geysircol[-1]) +
    ylab("Fusion calls (%)") +
    labs(fill = "Excluded by")
  plot(p)
}
#-----------------------------------

# TPM-FTS 2D Plot
create_TPM_FTS_ggplot <- function(rt,max_x,max_y,flabels){
  require(ggplot2)
  require(ggrepel)
  require(dplyr)
  
  rt_k <- rt[rt$known=="known",]
  rt_uk <- rt[rt$known=="unknown",]
  rt_uk <- rt_uk[sample(1:nrow(rt_uk),1000),]
  rt <- rbind.data.frame(rt_uk,rt_k)
  rt$filter_status <- factor("failed",levels=c("passed","failed"))
  rt$filter_status[rt$ev_level>3] <- "passed"
  rt <- arrange(rt,desc(known),FTS,filter_status)
  rt$known <- factor(rt$known,levels=c("known","unknown"))
  
  rt$tpm_gene1 <- log2(rt$tpm_gene1+1) 
  rt$tpm_gene2 <- log2(rt$tpm_gene2+1)
  
  #fusion coverage to total gene coverage ratio plot (length normalized)
  p <- ggplot(rt,aes(x=tpm_gene2,y=tpm_gene1)) +
    geom_point(aes(shape=known,fill=filter_status,size=FTS),alpha=.5) +
    scale_shape_manual(values=c(24,21)) +
    scale_size() +
    scale_fill_manual(values=c("#85D624","#FAA100")) +
    xlab("log2(TPM) 3' gene") + ylab("log2(TPM) 5' gene") +
    scale_x_continuous(breaks = pretty(0:max_x,n=10),limits = c(0,max_x)) +
    scale_y_continuous(breaks = pretty(0:max_y,n=10),limits = c(0,max_y)) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(),
          legend.justification = c("right","top"),
          legend.position = c(0.99,0.98),
          legend.box.just = "right",
          legend.spacing = unit(2,"char"),
          legend.background = element_rect(fill="white",color = "black", size = 0.2),
          legend.text = element_text(size=10),
          legend.title = element_text(size=11,face="bold")
    ) +
    guides(shape=guide_legend(order = 1, title = "Fusion event",
                              override.aes = list(size = 4)),
           size=guide_legend(order = 2, title = "FTS"),
           fill=guide_legend(order = 3, title = "FTS filter",
                             override.aes = list(shape=22,color="white",alpha=.8,
                                                 size = 5))
    )
  if(length(flabels)>0){
    rt$label <- ""
    rt[rt$fusion_name %in% flabels,]$label <- paste0(rt[rt$fusion_name %in% flabels,]$gene1,"-",rt[rt$fusion_name %in% flabels,]$gene2)
    p + geom_text_repel(data=rt,aes(label=label),force=15,size=1.5,segment.size = 0.2,
                        min.segment.length = 0.1,)
  }
  plot(p)
}
#-----------------------------------

# TPM-FTS 3D plot with plotly
create_TPM_FTS_plotly <- function(rt){
  require(plotly)
  require(dplyr)
  
  rt_k <- rt[rt$known=="known" & !rt$reciprocal,]
  rt_uk <- rt[rt$known=="unknown",]
  rt_uk <- rt_uk[sample(1:nrow(rt_uk),min(nrow(rt_uk),1000)),]
  rt <- rbind.data.frame(rt_uk,rt_k)
  rt$siz <- 1
  rt$siz[rt$known=="known"] <- 2
  
  rt$known <- factor(rt$known,levels=c("unknown","known"))
  rt$tpm5prime <- log2(rt$tpm5prime+1) 
  rt$tpm3prime <- log2(rt$tpm3prime+1)
  rt$tpmfusion <- log2(rt$tpmfusion+1)
  
  plot_ly(rt, x = ~tpm5prime, y = ~tpm3prime, z = ~FTS, 
          color = ~known, size = ~siz, alpha=1) %>%
    add_markers(colors = c("#ef8a62","#67a9cf"),
                sizes = c(70,210), symbols = c("circle-open"),
                stroke = I("black"), span = I(1),
                marker = list(opacity=1),
                hoverinfo = "text",
                text = ~paste0(" <b>Sample:</b>\t", sample,"<br>",
                               " <b>Fusion name:</b>\t", label,"<br>",
                               " <b>FTS:</b>\t\t", FTS,"<br><br>",
                               " <b>Expression in log2(TPM):</b><br>",
                               " ",gene1,":\t", tpm5prime,"<br>",
                               " ",gene2,":\t", tpm3prime,"<br>",
                               " Fusion transcript:\t",tpmfusion)) %>%
    layout(scene=list(xaxis=list(title='5\' Gene log2(TPM)'),
                      yaxis=list(title='3\' Gene log2(TPM)')),
           legend=list(title=list(text='<b> Fusion event </b>')))
}
#-----------------------------------

# Oncoprint
oncoprint <- function(outfile){
  if(all(is.na(clintable$Karyotype) | clintable$Karyotype == "") 
     | all(is.na(clintable$otherCyto) | clintable$otherCyto == "")){
    cat("\n    ",yellow("->")," No karyotype / molecular diagnostics information in clinical table. Oncoprint will be skipped.",sep="")
    return()
  }
  
  # Reduce mitelman table to 'known' fusion genes by means of ChimerDB
  # If ISCN karyo of fusion gene matches other fusion genes in database
  # Reduce Mitelman table to fusion genes with the highest recurrence
  mitelman_reduced <- mitelman[mitelman$fusion %in% chimerdb_known,]
  idx <- apply(mitelman_reduced, 1, function(x){
    mitelman_redundance_idx <- which(grepl(x['karyogrep'],mitelman_reduced$karyo))
    temp_idx <- which.max(mitelman_reduced$rec[mitelman_redundance_idx])
    mitelman_redundance_idx[temp_idx]
  })
  mitelman_reduced <- mitelman_reduced[idx,]
  
  # Check karyotypes in clinical table for occurence of aberrations reported to result in known fusion genes
  idx <- apply(mitelman_reduced, 1, function(x){
    temp_idx <- which(grepl(x['karyogrep'], clintable$Karyotype, ignore.case = T))
    if(length(temp_idx) > 0){
      data.frame(sample=clintable$sample[temp_idx], 
                 fusion=rep(x['fusion'],length(temp_idx)), 
                 cohort=rep(clintable$cohort[temp_idx],length(temp_idx)))
    }
  })
  karyo_search <- data.frame(sample=character(),fusion=character(),cohort=character(),karyo=numeric(),stringsAsFactors = F)
  if(!is.null(idx)){
    karyo_search <- do.call(rbind,idx)
    karyo_search <- unique(karyo_search)
    karyo_search$karyo <- 6  
  }

  # Check MDx in clinical table for reported known fusion genes
  idx <- apply(mitelman_reduced, 1, function(x){
    temp_idx <- which(grepl(x['fusion'], clintable$otherCyto, ignore.case = T))
    if(length(temp_idx) > 0){
      data.frame(sample=clintable$sample[temp_idx], 
                 fusion=rep(x['fusion'],length(temp_idx)), 
                 cohort=rep(clintable$cohort[temp_idx],length(temp_idx)))
    }
  })
  mol_search <- data.frame(sample=character(),fusion=character(),cohort=character(),mol=numeric(),stringsAsFactors = F)
  if(!is.null(idx)){
    mol_search <- do.call(rbind,idx)
    mol_search <- unique(mol_search)
    mol_search$mol <- 6
  }

  # Check RNA-seq for reported known fusion genes
  idx <- apply(mitelman_reduced, 1, function(x){
    temp_idx <- which(grepl(x['fusion'], resultTable$label, ignore.case = T))
    temp_idx <- c(temp_idx,which(grepl(x['fusion'], resultTable$reciproc_label, ignore.case = T)))
    if(length(temp_idx) > 0){
      data.frame(sample=resultTable$sample[temp_idx], 
                 fusion=rep(x['fusion'],length(temp_idx)), 
                 cohort=rep(resultTable$cohort[temp_idx],length(temp_idx)), 
                 rna=resultTable$ev_level[temp_idx])
    }
  })
  rna_search <- data.frame(sample=character(),fusion=character(),cohort=character(),rna=numeric(),stringsAsFactors = F)
  if(!is.null(idx)){
    rna_search <- do.call(rbind,idx)
    rna_search <- rna_search[order(rna_search$rna,decreasing = T),]
    rna_search <- rna_search[!duplicated(rna_search[,c("sample","fusion","cohort")]),] 
  }
  
  # Pre-processing of onco table for Oncoprint
  onco_matrix_long <- merge(karyo_search, mol_search, by=c("sample","fusion","cohort"),all = T)
  onco_matrix_long <- merge(onco_matrix_long, rna_search, by=c("sample","fusion","cohort"),all = T)
  onco_matrix_long[is.na(onco_matrix_long)] <- 0
  onco_matrix_long <- subset(onco_matrix_long,karyo==6 | mol==6)
  
  # Sorting table for oncoprint
  onco_matrix_long$score <- onco_matrix_long$karyo + onco_matrix_long$mol + onco_matrix_long$rna
  onco_matrix_long  <- onco_matrix_long %>%
    dplyr::group_by(fusion) %>%
    dplyr::mutate(rec=n()) 
  temp <- onco_matrix_long %>%
    dplyr::arrange(desc(rec),fusion,desc(score),desc(karyo),desc(mol),desc(rna),cohort) %>%
    dplyr::select(sample,fusion)
  col_order <- unique(temp$sample)
  
  # Define order for rows
  fus_anno <- unique(temp$fusion)
  row_order <- unlist(lapply(fus_anno,FUN = function(x){
    r <- c()
    r <- c(r,paste0("karyo_",x))
    r <- c(r,paste0("mol_",x))
    r <- c(r,paste0("rna_",x))
    return(r)
  }
  ))
  
  # Create onco matrix for ComplexHeatMap
  temp <- data.table::dcast(setDT(onco_matrix_long),
                            formula = sample~fusion, value.var = names(onco_matrix_long)[4:6])
  temp <- data.frame(temp,check.names = F)
  rownames(temp) <- temp$sample
  temp <- temp[,-1]
  onco_matrix_wide <- t(temp)  
  
  # Cohort labels
  cohort_labels <- onco_matrix_long$cohort[match(colnames(onco_matrix_wide),onco_matrix_long$sample)]
  
  # Row labels
  row_labels <- rownames(onco_matrix_wide)
  row_labels[grepl("karyo",row_labels)] <- "Karyotyping"
  row_labels[grepl("mol",row_labels)] <- "MDx"
  row_labels[grepl("rna",row_labels)] <- "RNA-seq"
  row_labels <- factor(row_labels,levels=c("Karyotyping","MDx","RNA-seq"))
  
  # Initializing Heatmap
  colors <- structure(c("grey","#edf8e9","#c7e9c0","#a1d99b","#74c476","#31a354","#006d2c"),
                      names=c("0","1","2","3","4","5","6"))
  
  num_cohorts <- length(unique(cohort_labels))
  cohort_cols <- RColorBrewer::brewer.pal(max(3,num_cohorts),"Set1")[1:num_cohorts]
  names(cohort_cols) <- unique(cohort_labels)
  
  text_param <- gpar(fontsize=18,col="black")
  title_param <- gpar(fontsize=18,col="black",fontface="bold")
  legend_param <- list(grid_height=unit(5,"mm"),grid_width=unit(5,"mm"),
                       title_gp=title_param,labels_gp=text_param)
  
  # Generate Heatmap
  p <- Heatmap(as.matrix(onco_matrix_wide),
               border = T,
               col = colors,
               rect_gp = gpar(col = "black", lwd = 0.5),
               na_col="grey",
               row_split = factor(gsub(".*_","",rownames(onco_matrix_wide)),levels=fus_anno),
               row_gap = unit(2, "mm"),
               row_title = NULL,
               cluster_row_slices = F,
               row_order = row_order,
               left_annotation = rowAnnotation(Fusion=anno_block(gp = gpar(fill=NULL,
                                                                           col=NA),
                                                                 labels = fus_anno,
                                                                 labels_gp = gpar(fontsize=16,col="black"),
                                                                 labels_rot = 0, labels_just = "right",
                                                                 labels_offset = unit(0.95, "npc"),
                                                                 width = unit(100, "mm")),
                                               Method=row_labels,
                                               annotation_legend_param = list(Method = c(legend_param,nrow=3)),
                                               col=list(Method=c("Karyotyping"="#283b42",
                                                                 "MDx"="#1d6a96",
                                                                 "RNA-seq"="#85b8cb")),
                                               show_annotation_name = F),
               column_order = col_order,
               show_row_names = F, column_names_gp = gpar(fontsize = 10),
               heatmap_legend_param = c(list(title = "Evidence",
                                             border = "black"),
                                        legend_param),
               top_annotation = HeatmapAnnotation(Cohort=cohort_labels,
                                                  annotation_legend_param = list(Cohort = c(legend_param,nrow=5)),
                                                  col=list(Cohort = cohort_cols),
                                                  show_annotation_name = F,
                                                  gp = gpar(col = "black")),
               width = ncol(onco_matrix_wide)*unit(4, "mm"), 
               height = nrow(onco_matrix_wide)*unit(4, "mm")
  )
  
  # Obtain dimensions of Heatmap    
  pdf(NULL)
  p <- draw(p,heatmap_legend_side = "right",annotation_legend_side = "top",
            adjust_annotation_extension = T, padding=unit(c(1,1,1,1),"mm")
            # heatmap_legend_list = list(Legend(labels = "low evidence", type = "points", pch = 8, background = NA))
  )
  w = ComplexHeatmap:::width(p)
  w = convertX(w, "points", valueOnly = TRUE)
  h = ComplexHeatmap:::height(p)
  h = convertY(h, "points", valueOnly = TRUE)
  dev.off()
  
  # Plot and save to file
  png(outfile, width = w+40, height = h+40)
  draw(p,heatmap_legend_side = "right",annotation_legend_side = "top",
      adjust_annotation_extension = T, padding=unit(c(1,1,1,1),"mm"))
  dev.off()
}
#-----------------------------------

# Circos plots
# Input is a table with columns:
# "label" = Fusion gene name
# "chr1","chr2","gene1","gene2"
# "known" = Known or Unknown fusion event. Known if found recurrently in ChimerDB and annotated as experimentally verified
# "rec" = Recurrence of the fusion in the cohort
# "mitelman_rec" = Recurrence of the fusion in the MitelmanDB
  circosPlotsCandidates <- function(rt){
    ch <- rt$cohort[1]
    status <- rt$known[1]
    
    rt <- distinct(rt,sample,label,.keep_all = T)
    
    links <- rt %>%
      select(label,break5prime,break3prime,mitelman_rec) %>%
      group_by(label) %>%
      mutate(rec=n())
    
    labels <- strsplit(links$label,"::")
    links$gene1 <- unlist(lapply(labels,function(x) x[1]))
    links$gene2 <- unlist(lapply(labels,function(x) x[2]))
    links <- distinct(links,label,.keep_all = T)

    temp <- do.call(rbind.data.frame, strsplit(links$break5prime,":"))
    links$chr1 <- gsub("chr","hs",temp[,1])
    links$start1 <- temp[,2]
    links$end1 <- links$start1
    temp <- do.call(rbind.data.frame, strsplit(links$break3prime,":"))
    links$chr2 <- gsub("chr","hs",temp[,1])
    links$start2 <- temp[,2]
    links$end2 <- links$start2
    
    links <- links[,match(c("chr1","start1","end1",
                            "chr2","start2","end2","rec","mitelman_rec","gene1","gene2","label"),colnames(links))]
    
    # Label all genes if status=="known", otherwise label only genes involved in recurrent fusions
      x <- 1
      if(status=="unknown") x <- 2
      labels <- rbind(setNames(links[links$rec>=x,match(c("chr1","start1","end1","gene1"),colnames(links))],c("chr","start","end","gene")),
                      setNames(links[links$rec>=x,match(c("chr2","start2","end2","gene2"),colnames(links))],c("chr","start","end","gene")))
      labels <- unique(labels)
    
    # Set line thickness and colors
      links$thickness <- rescale(as.numeric(links$rec),to = c(3,30))
      links$meta <- ""
      for(i in 1:nrow(links)){
        thickness <- round(links$thickness[i],digits = 0)
        recurrence <- links$rec[i]
        mitelman_recurrence <- links$mitelman_rec[i]
        links$meta[i] <- paste0("thickness=",thickness)
        if(status=="known"){
          links$meta[i] <- paste0(links$meta[i],",color=blue_a4,z=15")
        }
        else{
          if(recurrence>1){
            if(mitelman_recurrence>0){
              links$meta[i] <- paste0(links$meta[i],",color=green_a4,z=15")
            }else{
              links$meta[i] <- paste0(links$meta[i],",color=red_a4,z=20")
            }
          } else{
            if(mitelman_recurrence>0){
              links$meta[i] <- paste0(links$meta[i],",color=green_a4,z=15")
            }
          }
        }
      }
      links <- links[,-match(c("thickness","rec","mitelman_rec","gene1","gene2","label"),colnames(links))]
    
    system(paste0("mkdir -p ",outputfolder,"/filter_results/tmp"))
    write.table(links,paste0(outputfolder,"/filter_results/tmp/links.txt"),col.names = F,row.names = F,sep="\t",quote = F)
    write.table(labels,paste0(outputfolder,"/filter_results/tmp/labels.txt"),col.names = F,row.names = F,sep="\t",quote = F)
    
    cmd <- "/tools/circos-0.69-9/bin/circos"
    cmd <- paste0(cmd," -conf /scripts/circosPlotsConfig/circos.conf -silent")
    cmd <- paste0(cmd," -dir ",resultfolder,"/plots -file Circosplot_",ch,"_",status,"_fusions")
    system(cmd)
    unlink(paste0(outputfolder,"/filter_results/tmp/"),recursive = T)
  }
#-----------------------------------
