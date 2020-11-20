library(mixOmics)
library(clusterProfiler)
library(RDAVIDWebService)

setwd("your directory")
RESULTSDF = read.table("resultsdf.txt", header = TRUE, sep = "\t")


SAMPLES=RESULTSDF[,150:207]
SAMPLES_EDL=RESULTSDF[,150:172]
SAMPLES_HRT=RESULTSDF[,173:184]
SAMPLES_SOL=RESULTSDF[,185:207]

TYP=c('EDL_Tumour',     'EDL_Tumour_vib',  'EDL_Control_vib',  'EDL_Control',                          'EDL_Control_vib',
      'EDL_Tumour',     'EDL_Control',     'EDL_Control_vib',  'EDL_Control',      'EDL_Tumour_vib',   'EDL_Tumour',
      'EDL_Control',    'EDL_Tumour',      'EDL_Control_vib',  'EDL_Tumour_vib',   'EDL_Tumour_vib',   'EDL_Control_vib',
      'EDL_Control',    'EDL_Tumour',      'EDL_Tumour_vib',   'EDL_Control_vib',  'EDL_Tumour',       'EDL_Control',
      'HRT_Tumour_vib', 'HRT_Control',     'HRT_Control_vib',  
      'HRT_Tumour',     'HRT_Control_vib', 'HRT_Tumour_vib',
      'HRT_Control',    'HRT_Tumour_vib',  'HRT_Control',
      'HRT_Tumour',     'HRT_Control_vib', 'HRT_Tumour',
      'SOL_Tumour',     'SOL_Tumour_vib',  'SOL_Control_vib',  'SOL_Control',                          'SOL_Control_vib',
      'SOL_Tumour',     'SOL_Control',     'SOL_Control_vib',  'SOL_Control',      'SOL_Tumour_vib',   'SOL_Tumour',
      'SOL_Control',    'SOL_Tumour',      'SOL_Control_vib',  'SOL_Tumour_vib',   'SOL_Tumour_vib',   'SOL_Control_vib',  
      'SOL_Control',    'SOL_Tumour',      'SOL_Tumour_vib',   'SOL_Control_vib',  'SOL_Tumour',       'SOL_Control')

TIS=c('EDL','EDL','EDL','EDL',      'EDL',
      'EDL','EDL','EDL','EDL','EDL','EDL',
      'EDL','EDL','EDL','EDL','EDL','EDL',
      'EDL','EDL','EDL','EDL','EDL','EDL',
      'HRT','HRT','HRT',
      'HRT','HRT','HRT',
      'HRT','HRT','HRT',
      'HRT','HRT','HRT',
      'SOL','SOL','SOL','SOL',      'SOL',
      'SOL','SOL','SOL','SOL','SOL','SOL',
      'SOL','SOL','SOL','SOL','SOL','SOL',
      'SOL','SOL','SOL','SOL','SOL','SOL')

GRP=c('Tumour',     'Tumour_vib',  'Control_vib',  'Control',                      'Control_vib',
      'Tumour',     'Control',     'Control_vib',  'Control',      'Tumour_vib',   'Tumour',
      'Control',    'Tumour',      'Control_vib',  'Tumour_vib',   'Tumour_vib',   'Control_vib',
      'Control',    'Tumour',      'Tumour_vib',   'Control_vib',  'Tumour',       'Control',
      'Tumour_vib', 'Control',     'Control_vib',  
      'Tumour',     'Control_vib', 'Tumour_vib',
      'Control',    'Tumour_vib',  'Control',
      'Tumour',     'Control_vib', 'Tumour',
      'Tumour',     'Tumour_vib',  'Control_vib',  'Control',                      'Control_vib',
      'Tumour',     'Control',     'Control_vib',  'Control',      'Tumour_vib',   'Tumour',
      'Control',    'Tumour',      'Control_vib',  'Tumour_vib',   'Tumour_vib',   'Control_vib',  
      'Control',    'Tumour',      'Tumour_vib',   'Control_vib',  'Tumour',       'Control')

GRP=c('T',     'T+V',  'C+V',  'C',                      'C+V',
      'T',     'C',     'C+V',  'C',      'T+V',   'T',
      'C',    'T',      'C+V',  'T+V',   'T+V',   'C+V',
      'C',    'T',      'T+V',   'C+V',  'T',       'C',
      'T+V', 'C',     'C+V',  
      'T',     'C+V', 'T+V',
      'C',    'T+V',  'C',
      'T',     'C+V', 'T',
      'T',     'T+V',  'C+V',  'C',                      'C+V',
      'T',     'C',     'C+V',  'C',      'T+V',   'T',
      'C',    'T',      'C+V',  'T+V',   'T+V',   'C+V',  
      'C',    'T',      'T+V',   'C+V',  'T',       'C')

VIB=c('',     'vib',  'vib',  '',             'vib',
      '',     '',     'vib',  '',     'vib',  '',
      '',     '',     'vib',  'vib',  'vib',  'vib',
      '',     '',     'vib',  'vib',  '',     '',
      'vib',  '',     'vib',  
      '',     'vib',  'vib',
      '',     'vib',  '',
      '',     'vib',  '',
      '',     'vib',  'vib',  '',             'vib',
      '',     '',     'vib',  '',     'vib',  '',
      '',     '',     'vib',  'vib',  'vib',  'vib',  
      '',     '',     'vib',  'vib',  '',     '')

## EDL

plsda.res=splsda(t(SAMPLES_EDL),GRP[1:23],ncomp=2,keepX=rep(50,2),scale=TRUE)

tiff('EDL_PLSDA.tiff',width = 600, height = 400)
plotIndiv(plsda.res,
          group=GRP[1:23],
          comp = c(1,2),
          ind.names=TIS[1:23],
          ellipse=TRUE,
          title=paste("sPLS-DA EDL"),
          legend =TRUE,
          size.title = 20,
          size.xlabel=20,
          size.ylabel=20,
          size.axis=20,
          size.legend=15,
          size.legend.title = 15)
dev.off()

EDL_COMP2=selectVar(plsda.res,comp=2)
SYM=match(EDL_COMP2$name,rownames(RESULTSDF))
GENES <- cbind(data.frame(RESULTSDF$SYMBOL[SYM]),data.frame(RESULTSDF$ENTREZID[SYM]))
write.table(GENES, "PLSDA_EDL_COMP2.txt", quote=FALSE, col.names=NA, sep="\t")

EDL_COMP1=selectVar(plsda.res,comp=1)
SYM=match(EDL_COMP1$name,rownames(RESULTSDF))
GENES <- cbind(data.frame(RESULTSDF$SYMBOL[SYM]),data.frame(RESULTSDF$ENTREZID[SYM]))
write.table(GENES, "PLSDA_EDL_COMP1.txt", quote=FALSE, col.names=NA, sep="\t")

# ggo=enrichDAVID(gene = EDL_COMP2$name,idType = "ENTREZ_GENE_ID",listType = "Gene",annotation = "GOTERM_BP_ALL", david.user="rogier.plas@wur.nl")
# summary(data.frame(ggo))
# barplot(ggo)
# cnetplot(ggo)
# plot(ggo)

# library(org.Mm.eg.db)
# 
# # select 50 genes
# my.genes <- keys(org.Mm.eg.db)[1:50]
# GENES <- EDL_COMP2$name
# 
# # retrieve GO annotation info (Note thate multiple GO categories are present per gene)
# anno.result <- AnnotationDbi:::select(org.Mm.eg.db, keys=my.genes, columns=c("ENTREZID","GO"),keytype="ENTREZID")
# anno.result <- AnnotationDbi:::select(org.Mm.eg.db, keys=GENES, columns=c("ENTREZID","GO"),keytype="ENTREZID")
# 
# library(GO.db)
# # extract a named vector of all terms
# goterms <- Term(GOTERM)
# 
# # match GO ID with GO name
# m <- match(anno.result$GO, names(goterms))
# 
# # combine info
# GO.annotation <- cbind(anno.result, goterms[m])
# 
# # optional: filter/include specific terms
# # BP	Biological Process
# # IEA	Inferred from Electronic Annotation
# # ND	No biological Data available
# 
# 
# GO.annotation <- GO.annotation[GO.annotation$ONTOLOGY == "BP", ] #include only Biological Process
# 
# 
# #Combined filtering [exclude all IEA and ND (and 'NA')]
# GO.annotation <- GO.annotation[GO.annotation$EVIDENCE != "IEA" & 
#                                  GO.annotation$EVIDENCE != "ND" &
#                                  !is.na(GO.annotation$EVIDENCE)
#                                , ] 




## HRT

plsda.res=splsda(t(SAMPLES_HRT),GRP[24:35],ncomp=2,keepX=rep(50,2),scale=TRUE)

tiff('HRT_PLSDA.tiff',width = 600, height = 400)
plotIndiv(plsda.res,
          group=GRP[24:35],
          comp = c(1,2),
          ind.names=TIS[24:35],
          ellipse=TRUE,
          title=paste("sPLS-DA HRT"),
          legend =TRUE,
          size.title = 20,
          size.xlabel=20,
          size.ylabel=20,
          size.axis=20,
          size.legend=15,
          size.legend.title = 15)
dev.off()

HRT_COMP2=selectVar(plsda.res,comp=2)
SYM=match(HRT_COMP2$name,rownames(RESULTSDF))
GENES <- cbind(data.frame(RESULTSDF$SYMBOL[SYM]),data.frame(RESULTSDF$ENTREZID[SYM]))
write.table(GENES, "PLSDA_HRT_COMP2.txt", quote=FALSE, col.names=NA, sep="\t")

HRT_COMP1=selectVar(plsda.res,comp=1)
SYM=match(HRT_COMP1$name,rownames(RESULTSDF))
GENES <- cbind(data.frame(RESULTSDF$SYMBOL[SYM]),data.frame(RESULTSDF$ENTREZID[SYM]))
write.table(GENES, "PLSDA_HRT_COMP1.txt", quote=FALSE, col.names=NA, sep="\t")

## SOL

plsda.res=splsda(t(SAMPLES_SOL),GRP[36:58],ncomp=2,keepX=rep(50,2),scale=TRUE)

tiff('SOL_PLSDA.tiff',width = 600, height = 400)
plotIndiv(plsda.res,
          group=GRP[36:58],
          comp = c(1,2),
          ind.names=TIS[36:58],
          ellipse=TRUE,
          title=paste("sPLS-DA SOL"),
          legend =TRUE,
          size.title = 20,
          size.xlabel=20,
          size.ylabel=20,
          size.axis=20,
          size.legend=15,
          size.legend.title = 15)
dev.off()

SOL_COMP2=selectVar(plsda.res,comp=2)
SYM=match(SOL_COMP2$name,rownames(RESULTSDF))
GENES <- cbind(data.frame(RESULTSDF$SYMBOL[SYM]),data.frame(RESULTSDF$ENTREZID[SYM]))
write.table(GENES, "PLSDA_SOL_COMP2.txt", quote=FALSE, col.names=NA, sep="\t")

SOL_COMP1=selectVar(plsda.res,comp=1)
SYM=match(SOL_COMP1$name,rownames(RESULTSDF))
GENES <- cbind(data.frame(RESULTSDF$SYMBOL[SYM]),data.frame(RESULTSDF$ENTREZID[SYM]))
write.table(GENES, "PLSDA_SOL_COMP1.txt", quote=FALSE, col.names=NA, sep="\t")

# DRAW VENN DIAGRAMS
COMP1_SOL_EDL_HRT=length(intersect(intersect(SOL_COMP1$name,EDL_COMP1$name),HRT_COMP1$name))
COMP1_SOL_HRT=length(intersect(SOL_COMP1$name,HRT_COMP1$name))
COMP1_SOL_EDL=length(intersect(SOL_COMP1$name,EDL_COMP1$name))
COMP1_EDL_HRT=length(intersect(EDL_COMP1$name,HRT_COMP1$name))
COMP1_EDL=length(EDL_COMP1$name)
COMP1_SOL=length(SOL_COMP1$name)
COMP1_HRT=length(EDL_COMP1$name)

library(VennDiagram)

grid.newpage()
tiff('COMP1_PLSDA.tiff',width = 500, height = 500)
VENN <- draw.triple.venn(area1 = COMP1_EDL, 
                 area2 = COMP1_SOL, 
                 area3 = COMP1_HRT, 
                 n12 = COMP1_SOL_EDL, 
                 n23 = COMP1_SOL_HRT, 
                 n13 = COMP1_EDL_HRT, 
                 n123 = COMP1_SOL_EDL_HRT, 
                 scaled = FALSE, 
                 euler.d = FALSE, 
                 category = c("EDL","SOL","HRT"), 
                 lty = "blank", 
                 fill = c("green","blue","yellow"), 
                 print.mode = as.character(c('raw','percent')),
                 fontfamily="Arial",
                 fontface="plain",
                 cat.fontfamily="Arial",
                 cat.cex=2,
                 cex=1.75)
dev.off()

COMP2_SOL_EDL_HRT=length(intersect(intersect(SOL_COMP2$name,EDL_COMP2$name),HRT_COMP2$name))
COMP2_SOL_HRT=length(intersect(SOL_COMP2$name,HRT_COMP2$name))
COMP2_SOL_EDL=length(intersect(SOL_COMP2$name,EDL_COMP2$name))
COMP2_EDL_HRT=length(intersect(EDL_COMP2$name,HRT_COMP2$name))
COMP2_EDL=length(EDL_COMP2$name)
COMP2_SOL=length(SOL_COMP2$name)
COMP2_HRT=length(EDL_COMP2$name)


tiff('COMP2_PLSDA.tiff',width = 500, height = 500)
VENN <- draw.triple.venn(area1 = COMP2_EDL, 
                         area2 = COMP2_SOL, 
                         area3 = COMP2_HRT, 
                         n12 = COMP2_SOL_EDL, 
                         n23 = COMP2_SOL_HRT, 
                         n13 = COMP2_EDL_HRT, 
                         n123 = COMP2_SOL_EDL_HRT, 
                         scaled = FALSE, 
                         euler.d = FALSE, 
                         category = c("EDL","SOL","HRT"), 
                         lty = "blank", 
                         fill = c("green","blue","yellow"), 
                         print.mode = as.character(c('raw','percent')),
                         fontfamily="Arial",
                         fontface="plain",
                         cat.fontfamily="Arial",
                         cat.cex=2,
                         cex=1.75)
dev.off()