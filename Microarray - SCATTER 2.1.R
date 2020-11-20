SCATTER <- function(FC,SIGN,X,Y,int,slop,Xname,Yname){
  
  # FC = fold change matrix
  # SIGN = TRUE/FALSE matrix with significances corresponding to fold change values of FC
  # X = colname of X variable in both FC and SIGN matrix
  # Y = colname of Y variable in both FC and SIGN matrix

  COLOR=ifelse((SIGN[,X]==TRUE & SIGN[,Y]==TRUE),"BOTH",ifelse((SIGN[,X]==TRUE),"X",ifelse((SIGN[,Y]==TRUE),"Y","NONE")))

  TEMP=cbind(FC,COLOR)
  
  plot = ggplot(FC, aes(x=eval(parse(text=colnames(FC[X]))), y=eval(parse(text=colnames(FC[Y]))), color=COLOR))
  
  plot +
    geom_point(size = 3) +
    theme_bw() +
    geom_vline(xintercept = 0) + # plot vertical line
    geom_hline(yintercept = 0) + # plot horizontal line
    geom_abline(intercept = 0, slope = 1,size=0.3,colour="black") +
    geom_abline(intercept = int, slope = slop, size = 1,colour="orange") +
    scale_color_manual(name="Significance: ",
                       breaks=c("NONE", "X", "Y", "BOTH"),
                       label = c("NONE",Xname,Yname,"BOTH"),
                       values = c('#377EB8','black', '#E41A1C', '#4DAF4A')) +
    xlab(Xname) + 
    ylab(Yname) + 
    coord_cartesian(ylim = c(-6,6), xlim = c(-6,6)) +
    # coord_fixed(ratio = 1) +
    theme(strip.text.x = element_text(size=40, face="bold"),
          axis.text=element_text(size=16),
          axis.title=element_text(size=16,face=c('bold')),
          axis.title.y = element_text(vjust =1.5),
          legend.title = element_text(size=16, face="bold.italic"),
          legend.text = element_text(size = 16),
          legend.key = element_blank(),
          legend.position="bottom")
}
