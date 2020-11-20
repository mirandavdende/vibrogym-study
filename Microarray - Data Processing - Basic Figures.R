library(RColorBrewer)
library(ggplot2)
library(extrafont)
library(VennDiagram)
font_import()
loadfonts(device = "win")

setwd("your directory")
RESULTSDF <- read.table("resultsdf.txt", header = TRUE, sep = "\t")


## SET PARAMS

ADJ.P=0.01
PHENODATA <- data.frame(cbind(rep(c("green3","darkblue"),3),rep(c("No Training","WBV"),3)))
colnames(PHENODATA) <- c("COLORCODE","TYPE")

## SELECT 

P=data.frame(cbind(RESULTSDF$adj.IBMT.p.EDLCvsT....EDL_T...EDL_C<ADJ.P,RESULTSDF$adj.IBMT.p.EDLCvsTV....EDL_TV...EDL_C<ADJ.P,RESULTSDF$adj.IBMT.p.HRTCvsT....HRT_T...HRT_C<ADJ.P,RESULTSDF$adj.IBMT.p.HRTCvsTV....HRT_TV...HRT_C<ADJ.P,RESULTSDF$adj.IBMT.p.SOLCvsT....SOL_T...SOL_C<ADJ.P,RESULTSDF$adj.IBMT.p.SOLCvsTV....SOL_TV...SOL_C<ADJ.P))
FC=data.frame(cbind(RESULTSDF$coefficients.EDLCvsT....EDL_T...EDL_C,RESULTSDF$coefficients.EDLCvsTV....EDL_TV...EDL_C,RESULTSDF$coefficients.HRTCvsT....HRT_T...HRT_C,RESULTSDF$coefficients.HRTCvsTV....HRT_TV...HRT_C,RESULTSDF$coefficients.SOLCvsT....SOL_T...SOL_C,RESULTSDF$coefficients.SOLCvsTV....SOL_TV...SOL_C))
colnames(P) <- c("EDL No training","EDL WBV","Heart No training","Heart WBV","SOL No training","SOL WBV")
colnames(FC) <- c("EDL No training","EDL WBV","Heart No training","Heart WBV","SOL No training","SOL WBV")

## VENN DIAGRAMS
source('M:/WUR/Projects/RO2.1 - Mouse Fingerprint (Experimental C26)/Measurements/150201 - 0305 Main Experiment/Microarray/4. R_codes/171124 VENNY 2.2.R')
setEPS()
postscript('EDL_VENN_T_vs_NOT.eps',width = 600, height = 600,family = "Arial")
VENNY(P,c(1,2),PHENODATA)
dev.off()
tiff('HRT_VENN_T_vs_NOT.tiff',width = 600, height = 600)
VENNY(P,c(3,4),PHENODATA)
dev.off()
tiff('SOL_VENN_T_vs_NOT.tiff',width = 600, height = 600)
VENNY(P,c(5,6),PHENODATA)
dev.off()

colnames(P) <- c("EDL_C_vs_T","EDL_C_vs_TV","HRT_C_vs_T","HRT_C_vs_TV","SOL_C_vs_T","SOL_C_vs_TV")
colnames(FC) <- c("EDL_C_vs_T","EDL_C_vs_TV","HRT_C_vs_T","HRT_C_vs_TV","SOL_C_vs_T","SOL_C_vs_TV")


## SCATTER PLOTS
source('M:/WUR/Projects/RO2.1 - Mouse Fingerprint (Experimental C26)/Measurements/150201 - 0305 Main Experiment/Microarray/4. R_codes/170822 SCATTER 2.1.R')
tiff('EDL_SCATTER_T_vs_NOT.tiff',width = 600, height = 600)
EDL_coef=lm(EDL_C_vs_TV ~ EDL_C_vs_T, data=FC)
SCATTER(FC,P,"EDL_C_vs_T","EDL_C_vs_TV",as.numeric(EDL_coef$coefficients[1]),as.numeric(EDL_coef$coefficients[2]),"EDL No training","EDL WBV")
cor(FC$EDL_C_vs_T,FC$EDL_C_vs_TV)
dev.off()
tiff('HRT_SCATTER_T_vs_NOT.tiff',width = 600, height = 600)
HRT_coef=lm(HRT_C_vs_TV ~ HRT_C_vs_T, data=FC)
SCATTER(FC,P,"HRT_C_vs_T","HRT_C_vs_TV",as.numeric(HRT_coef$coefficients[1]),as.numeric(HRT_coef$coefficients[2]),"Heart No training","Heart WBV")
cor(FC$HRT_C_vs_T,FC$HRT_C_vs_TV)
dev.off()
tiff('SOL_SCATTER_T_vs_NOT.tiff',width = 600, height = 600)
SOL_coef=lm(SOL_C_vs_TV ~ SOL_C_vs_T, data=FC)
SCATTER(FC,P,"SOL_C_vs_T","SOL_C_vs_TV",as.numeric(SOL_coef$coefficients[1]),as.numeric(SOL_coef$coefficients[2]),"SOL No training","SOL WBV")
cor(FC$SOL_C_vs_T,FC$SOL_C_vs_TV)
dev.off()

