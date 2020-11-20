library(stringr)
library(ggplot2)
library(forcats)
library(clusterProfiler)

setwd("your directory")
load("GSEA.Rdata")



names(GSEA)[c(1,4,6,9,11,14)] <-c("EDL C vs T","EDL C vs T+V","SOL C vs T","SOL C vs T+V","HRT C vs T","HRT C vs T+V")

merged.res <- as.data.frame(merge_result(GSEA[c(1,4,6,9,11,14)]))


## Select genesets with p<p_co
p_co=0.005
merged.res <- subset(merged.res, p.adjust<p_co)

## Select genesets which occur more that xxx
xxx=0
tt=table(merged.res$ID)
merged.res <- subset(merged.res, ID %in% names(tt[tt > xxx]))

## Remove genesets with disease types 
diseases=c("Myocarditis","leukemia","Measles","Leishmaniasis","disease","Pertussis","Malaria","infection","Tuberculosis","Legionellosis","Amphetamine")
merged.res <- merged.res[!grepl(paste(diseases,collapse="|"),merged.res$Description),]

## Set up/downregulation
merged.res$type = "upregulated"
merged.res$type[merged.res$NES < 0] = "downregulated"
merged.res$muscle = substr(merged.res$Cluster,1,3)
merged.res$muscle=factor(merged.res$muscle,levels = c("EDL","SOL","HRT"))
merged.res$vib = substr(merged.res$Cluster,11,12)
VIB=as.factor(c(merged.res$vib,"No Training","Training"))
VIB[VIB == ""] <- "No Training";VIB[VIB == "+V"] <- "Training";
merged.res$vib = VIB[1:107]

# Add Count and GeneRatio info
merged.res$Count <- str_count(merged.res$core_enrichment, "/") + 1
merged.res$GeneRatio <- merged.res$Count / merged.res$setSize

# Plot
p <- ggplot(merged.res, aes(x = NES, y = fct_reorder(Description , NES))) + 
  geom_point(aes(size = setSize, color = pvalue ,shape = factor(Cluster))) +
  # scale_colour_gradient(limits=c(0, p_co+0.001), low="red", breaks=c(0,p_co, 0.01),labels=c(0,p_co,0.01)) +
  theme_bw(base_size = 10) +
  xlab("Normalized Enrichment Score") +
  ylab(NULL) +
  ggtitle("KEGG pathway enrichment")

p

# Plot
p <- ggplot(merged.res, aes(x = NES, y = fct_reorder(Description , NES))) + 
  geom_point(aes(size = setSize, color = p.adjust, shape = factor(vib))) +
  scale_shape_manual(values=c(16,2)) +
  scale_color_gradient(limits=c(0, p_co+0.001), low="red", breaks=c(0,p_co, 0.01),labels=c(0,p_co,0.01)) +
  theme_bw(base_size = 10) +
  xlab("Normalized Enrichment Score") +
  ylab(NULL) +
  ggtitle("KEGG pathway enrichment")

p

# For PDF
pp <- p + facet_grid(.~muscle) #Cluster and type are columns to 'split' on
jpeg('GSEA.jpeg',width = 1000, height = 600,quality = 100)
pp
dev.off()


merged.res.up=subset(merged.res, merged.res$type == "upregulated")
merged.res.down=subset(merged.res, merged.res$type == "downregulated")

# Plot
p <- ggplot(merged.res.up, aes(x = NES, y = fct_reorder(Description , NES))) + 
  geom_point(aes(size = setSize, color = p.adjust)) +
  scale_colour_gradient(limits=c(0, p_co+0.001), low="red", breaks=c(0,p_co, 0.01),labels=c(0,p_co,0.01)) +
  theme_bw(base_size = 10) +
  xlab("Normalized Enrichment Score") +
  ylab(NULL) +
  ggtitle("KEGG pathway enrichment - Upregulation")

p

# For PDF
pp <- p + facet_grid(.~Cluster) #Cluster and type are columns to 'split' on
pp

# Plot
p <- ggplot(merged.res.down, aes(x = NES, y = fct_reorder(Description , -NES))) + 
  geom_point(aes(size = setSize, color = p.adjust)) +
  scale_colour_gradient(limits=c(0, p_co+0.001), low="red", breaks=c(0,p_co, 0.01),labels=c(0,p_co,0.01)) +
  theme_bw(base_size = 10) +
  xlab("Normalized Enrichment Score") +
  ylab(NULL) +
  ggtitle("KEGG pathway enrichment - Downregulation")

p

# For PDF
pp <- p + facet_grid(.~Cluster) #Cluster and type are columns to 'split' on
pp
# ggsave(filename="Absolute_GSEA.NES.KEGG.pdf", plot=pp, width = 21, height = 29.1, units = "cm")
