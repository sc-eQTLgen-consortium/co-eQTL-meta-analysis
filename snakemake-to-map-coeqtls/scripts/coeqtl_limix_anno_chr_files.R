

## Generate files per chromosome

ct = 'Mono'

annot <- read.delim(paste0("./LimixAnnotationFile_co_",ct,".txt.gz"))
e_annot <- read.delim(paste0("./LimixAnnotationFile_co2_",ct,".txt.gz"))
fvf <-  read.delim(paste0("../all_features_to_test_",ct,".tsv.gz"))

for (chr in 22:1){
  print(paste("chromsome",chr))
  annotRel <- annot[which(annot$chromosome==chr),]
  e_annotRel <- e_annot[which(e_annot$chromosome==chr),]
  
  annotRel <- rbind(annotRel,e_annotRel)
  
  fvfRel = fvf[which(fvf$feature_id %in% annotRel$feature_id),]
  
  e_annotRel = annotRel[which(duplicated(annotRel$feature_id)),]
  annotRel = annotRel[which(!duplicated(annotRel$feature_id)),]
  
  write.table(annotRel,file = paste0("./LimixAnnotationFile_co_",ct,"_chr",chr,".txt"),sep="\t",quote=F,row.names=F)
  write.table(e_annotRel,file = paste0("./LimixAnnotationFile_co2_",ct,"_chr",chr,".txt"),sep="\t",quote=F,row.names=F)
  write.table(fvfRel,file = paste0("./all_features_to_test_",ct,"_chr",chr,".txt"),sep="\t",quote=F,row.names=F)
}


