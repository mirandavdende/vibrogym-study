# GSEA

setwd("your directory")
RESULTSDF = read.table("resultsdf.txt", header = TRUE, sep = "\t")

library(clusterProfiler)

GSEA=list()
IBMT.t_cols<-colnames(RESULTSDF)[grep("^IBMT.t.",colnames(RESULTSDF))]

for(i in 1:length(IBMT.t_cols)){
  genelist <- RESULTSDF[,IBMT.t_cols[i]];
  names(genelist)<-RESULTSDF$ENTREZID;
  genelist<-sort(genelist,decreasing = TRUE);
  GSEA[[i]] <- gseKEGG(geneList  = genelist,organism = 'mmu',nPerm = 20000, minGSSize = 10,pvalueCutoff = 1,verbose = FALSE)
}

names(GSEA)<-substr(IBMT.t_cols,8,nchar(IBMT.t_cols));
save(GSEA, file="GSEA.RData")


