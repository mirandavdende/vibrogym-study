VENNY <- function(SIG,SELECT,PHENODATA) {
    
  # COLOR=c("blue", "red", "yellow","green") #("skyblue", "pink1", "mediumorchid")
  COLOR=as.character(PHENODATA$COLORCODE[SELECT]);
  TYP=PHENODATA$TYPE[SELECT];  
  TIS=PHENODATA$TISSUE[SELECT];
  EXP=PHENODATA$ExperimentNr[SELECT];
  names=paste(TYP,"\n",EXP)
  
  # print(paste(Experiments, collapse = "|"));
  # print(grepl(paste(Experiments, collapse = "|"),PHENODATA$Experiment));
  # print(names);
  # print(Experiments);
  grid.newpage()
  if (length(SELECT) == 1) {
    draw.single.venn(area = nrow(subset(SIG, SIG[,SELECT[1]]== TRUE)),
                            category = names, 
                            lty = "blank", 
                            fill = COLOR[1],
                            theme(text=element_text(size=16,family="TT Arial")))
    
  }
  
  if (length(SELECT) == 2) {
    draw.pairwise.venn(area1 = nrow(subset(SIG, SIG[,SELECT[1]]== TRUE)), 
                            area2 = nrow(subset(SIG, SIG[,SELECT[2]]== TRUE)), 
                            cross.area = nrow(subset(SIG, SIG[,SELECT[1]]== TRUE & SIG[,SELECT[2]]== TRUE)), 
                            category = names, 
                            lty = "blank", 
                            fill = COLOR[1:2],
                            scaled = FALSE, 
                            print.mode = as.character(c('raw','percent')),
                            fontfamily="Arial",
                            fontface="plain",
                            cat.fontfamily="Arial",
                            cat.cex=2,cex=2,
                            cat.pos=1)
  }
  
  if (length(SELECT) == 3) {
    draw.triple.venn(area1 = nrow(subset(SIG, SIG[,SELECT[1]]== TRUE)), 
                            area2 = nrow(subset(SIG, SIG[,SELECT[2]]== TRUE)), 
                            area3 = nrow(subset(SIG, SIG[,SELECT[3]]== TRUE)), 
                            n12 = nrow(subset(SIG, SIG[,SELECT[1]]== TRUE & SIG[,SELECT[2]]== TRUE)), 
                            n23 = nrow(subset(SIG, SIG[,SELECT[2]]== TRUE & SIG[,SELECT[3]]== TRUE)), 
                            n13 = nrow(subset(SIG, SIG[,SELECT[1]]== TRUE & SIG[,SELECT[3]]== TRUE)), 
                            n123 = nrow(subset(SIG, SIG[,SELECT[1]]== TRUE & SIG[,SELECT[2]]== TRUE & SIG[,SELECT[3]]== TRUE)), 
                            category = names, 
                            lty = "blank", 
                            fill = COLOR[1:3])
  }
  
  if (length(SELECT) == 4) {
    draw.quad.venn(area1 = nrow(subset(SIG, SIG[,SELECT[1]]== TRUE)), 
                          area2 = nrow(subset(SIG, SIG[,SELECT[2]]== TRUE)), 
                          area3 = nrow(subset(SIG, SIG[,SELECT[3]]== TRUE)), 
                          area4 = nrow(subset(SIG, SIG[,SELECT[4]]== TRUE)), 
                          n12 = nrow(subset(SIG, SIG[,SELECT[1]]== TRUE & SIG[,SELECT[2]]== TRUE)), 
                          n23 = nrow(subset(SIG, SIG[,SELECT[2]]== TRUE & SIG[,SELECT[3]]== TRUE)), 
                          n34 = nrow(subset(SIG, SIG[,SELECT[3]]== TRUE & SIG[,SELECT[4]]== TRUE)),
                          n13 = nrow(subset(SIG, SIG[,SELECT[1]]== TRUE & SIG[,SELECT[3]]== TRUE)), 
                          n14 = nrow(subset(SIG, SIG[,SELECT[1]]== TRUE & SIG[,SELECT[4]]== TRUE)), 
                          n24 = nrow(subset(SIG, SIG[,SELECT[2]]== TRUE & SIG[,SELECT[4]]== TRUE)), 
                          n123 = nrow(subset(SIG, SIG[,SELECT[1]]== TRUE & SIG[,SELECT[2]]== TRUE & SIG[,SELECT[3]]== TRUE)), 
                          n234 = nrow(subset(SIG, SIG[,SELECT[2]]== TRUE & SIG[,SELECT[3]]== TRUE & SIG[,SELECT[4]]== TRUE)), 
                          n134 = nrow(subset(SIG, SIG[,SELECT[1]]== TRUE & SIG[,SELECT[3]]== TRUE & SIG[,SELECT[4]]== TRUE)), 
                          n124 = nrow(subset(SIG, SIG[,SELECT[1]]== TRUE & SIG[,SELECT[2]]== TRUE & SIG[,SELECT[4]]== TRUE)), 
                          n1234 = nrow(subset(SIG, SIG[,SELECT[1]]== TRUE & SIG[,SELECT[2]]== TRUE & SIG[,SELECT[3]]== TRUE & SIG[,SELECT[4]]== TRUE)),
                          category = names,
                          lty = "blank", 
                          fill = COLOR[1:4])
    
    
  }
}
