# Check whether all files are available for the samples listed in clintable
  file_missing = F
  for(s in clintable$sample){
    # check arriba calls
      if(!file.exists(paste0(outputfolder,"/arriba/",s,"/fusions.tsv"))){
        file_missing <- T
        cat(red("ERROR: "),"Missing file \".../arriba/",s,"/fusions.tsv\"\n",sep="")
        break
      }
    # check fusioncatcher calls
      if(!file.exists(paste0(outputfolder,"/fusioncatcher/",s,"/final-list_candidate-fusion-genes.txt"))){
        file_missing <- T
        cat(red("ERROR: "),"Missing file \".../fusioncatcher/",s,"/final-list_candidate-fusion-genes.txt\"\n",sep="")
        break
      }
    # check featurecounts
      if(!file.exists(paste0(outputfolder,"/featurecounts/reformatted/",s,".fc"))
         || file.size(paste0(outputfolder,"/featurecounts/reformatted/",s,".fc")) < 1000){
        file_missing <- T
        cat(red("ERROR: "),"Missing or truncated file \".../featurecounts/reformatted/",s,".fc\"\n",sep="")
        break
      }
    # check insertsizes
      if(!file.exists(paste0(outputfolder,"/insertsizes/",s,"_is_metrics.txt"))
         || file.size(paste0(outputfolder,"/insertsizes/",s,"_is_metrics.txt")) < 1000){
        file_missing <- T
        cat(red("ERROR: "),"Missing or truncated file \".../insertsizes/",s,"_is_metrics.txt\"\n",sep="")
        break
      }    
  }
  
  if(file_missing){
    print("Exiting analysis.")
    quit(save="no",status=0,runLast=F)
  }