# This is the R script used to analyse the data of the vibrogym study. 
# Shapiro Wilk, skewness and Kurtosis analysis was done on the whole dataset.
# Next, a Random Forest analysis and PCA analysis were preformed.
# R version 4.0.2 was used

# Install and load Packages
install.packages("moments")
library("moments") # Used for skewness and kurtosis

install.packages("ggplot2") # For making figures
library("ggplot2")

if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/ggpubr")
install.packages("ggpubr") # For making figure layouts
library("ggpubr")

install.packages("randomForest") # Used for random forest analysis
library("randomForest")

# Set working directory
setwd("your directory")

# Load data file
DF <- read.table("dataset.txt", header = TRUE)

# Transpose the data
DF.t <- t(DF)
DF.t2 <- as.data.frame(DF.t)

# Shapiro Wilk, skewness and Kurtosis analysis on the whole dataset.
lshap <- lapply(DF.t2, shapiro.test)
lskew <- lapply(DF.t2, skewness)
lkurt <- lapply(DF.t2, kurtosis)

# Extract the statistic and p.value
lres <- sapply(lshap, `[`, c("statistic","p.value"))
Temp.Table1 <- data.frame(t(lres))

# Make table including YES or NO for normality
SP.Table <- data.frame ("ShapiroWilk.p.value" = 0, "Normal" = 0)
SP.Table[nrow(SP.Table)+62,] <- NA
rownames(SP.Table) <- c(colnames(DF.t2))
SP.Table$ShapiroWilk.p.value <- Temp.Table1$p.value
SP.Table[,"Normal"] <- 'YES'
for (row in 1:nrow(SP.Table)) {
  if(SP.Table[row,"ShapiroWilk.p.value"] < 0.05) {SP.Table[row,"Normal"] <- 'NO'}
}
SP.Table$Skewness <- lskew
SP.Table$Kurtosis <- lkurt
SP.Table <- format(SP.Table, digits = 2)
SP.Table

# Check NA samples in the dataframe
sum(is.na(DF.t2))
which(!complete.cases(DF.t2))
colSums(is.na(DF.t2))

# Replace all NA values with mean
DF.t3 <- DF.t2
#Replace 1 value of Brains with mean
DF.t3$Brains [which(is.na(DF.t2$Brains ))] <- mean(DF.t2$Brains, na.rm = TRUE) 
# Replace 1 value of Atrogin.mRNA.exp with mean
DF.t3$Atrogin.mRNA.exp[which(is.na(DF.t2$Atrogin.mRNA))] <- mean(DF.t2$Atrogin.mRNA.exp, na.rm = TRUE)
# Replace 2 values of Daily.Activity with mean
DF.t3$Daily.Activity[which(is.na(DF.t2$Daily.Activity))] <- mean(DF.t2$Daily.Activity, na.rm = TRUE)
# Replace 2 values of Plasma.IFNgamma with mean
DF.t3$Plasma.IFNgamma [which(is.na(DF.t2$Plasma.IFNgamma ))] <- mean(DF.t2$Plasma.IFNgamma, na.rm = TRUE)
# Replace 2 values of Plasma.IL1beta with mean
DF.t3$Plasma.IL1beta [which(is.na(DF.t2$Plasma.IL1beta ))] <- mean(DF.t2$Plasma.IL1beta, na.rm = TRUE)
# Replace 2 values of Plasma.IL4 with mean
DF.t3$Plasma.IL4 [which(is.na(DF.t2$Plasma.IL4 ))] <- mean(DF.t2$Plasma.IL4, na.rm = TRUE)
# Replace 2 values of Plasma.IL6 with mean
DF.t3$Plasma.IL6 [which(is.na(DF.t2$Plasma.IL6 ))] <- mean(DF.t2$Plasma.IL6, na.rm = TRUE)
# Replace 2 values of Plasma.IL10 with mean
DF.t3$Plasma.IL10 [which(is.na(DF.t2$Plasma.IL10 ))] <- mean(DF.t2$Plasma.IL10, na.rm = TRUE)
# Replace 2 values of Plasma.MCP1 with mean
DF.t3$Plasma.MCP1 [which(is.na(DF.t2$Plasma.MCP1 ))] <- mean(DF.t2$Plasma.MCP1, na.rm = TRUE)
# Replace 2 values of Plasma.TNFalpha with mean
DF.t3$Plasma.TNFalpha [which(is.na(DF.t2$Plasma.TNFalpha ))] <- mean(DF.t2$Plasma.TNFalpha, na.rm = TRUE)
# Replace 2 values of Plasma.VEGF with mean
DF.t3$Plasma.VEGF [which(is.na(DF.t2$Plasma.VEGF ))] <- mean(DF.t2$Plasma.VEGF, na.rm = TRUE)
# Replace 1 value of CSA.SOL with mean
DF.t3$CSA.SOL [which(is.na(DF.t2$CSA.SOL ))] <- mean(DF.t2$CSA.SOL, na.rm = TRUE)
# Replace 1 value of CSA.fiber.SOL with mean
DF.t3$CSA.fiber.SOL [which(is.na(DF.t2$CSA.fiber.SOL ))] <- mean(DF.t2$CSA.fiber.SOL, na.rm = TRUE)
# Replace 1 value of CSA.typeI.SOL with mean
DF.t3$CSA.typeI.SOL [which(is.na(DF.t2$CSA.typeI.SOL ))] <- mean(DF.t2$CSA.typeI.SOL, na.rm = TRUE)
# Replace 1 value of CSA.typeII.SOL with mean
DF.t3$CSA.typeII.SOL [which(is.na(DF.t2$CSA.typeII.SOL ))] <- mean(DF.t2$CSA.typeII.SOL, na.rm = TRUE)
# Replace 1 value of CSA.ratio.SOL with mean
DF.t3$CSA.ratio.SOL [which(is.na(DF.t2$CSA.ratio.SOL ))] <- mean(DF.t2$CSA.ratio.SOL, na.rm = TRUE)
# Replace 1 value of AbAb.typeI.SOL with mean
DF.t3$AbAb.typeI.SOL [which(is.na(DF.t2$AbAb.typeI.SOL ))] <- mean(DF.t2$AbAb.typeI.SOL, na.rm = TRUE)
# Replace 1 value of RelAb.typeI.SOL with mean
DF.t3$RelAb.typeI.SOL [which(is.na(DF.t2$RelAb.typeI.SOL ))] <- mean(DF.t2$RelAb.typeI.SOL, na.rm = TRUE)
# Replace 1 value of AbAb.typeII.SOL with mean
DF.t3$AbAb.typeII.SOL [which(is.na(DF.t2$AbAb.typeII.SOL ))] <- mean(DF.t2$AbAb.typeII.SOL, na.rm = TRUE)
# Replace 1 value of RelAb.typeII.SOL with mean
DF.t3$RelAb.typeII.SOL [which(is.na(DF.t2$RelAb.typeII.SOL ))] <- mean(DF.t2$RelAb.typeII.SOL, na.rm = TRUE)

# Check if all the NA values are gone
sum(is.na(DF.t3))

#########################################################################################################

# Random Forest analysis
# Set seed so we can reproduce the results
set.seed(42)

# Create a new dataframe to add the groups to the data
Group <- data.frame("Group" = c("C","C","C","C","C","C","T","T","T","T","T","T",
                                "T+WBV","T+WBV","T+WBV","T+WBV","T+WBV","C+WBV","C+WBV","C+WBV","C+WBV","C+WBV","C+WBV"))

# Combine groups to data frame
DF.groups <- cbind(DF.t3, Group)

# Make the random forest model
model <- randomForest(as.factor(Group) ~ ., data = DF.groups, ntree=1000,  mtry =8, proximity=TRUE)
model

# Plot the error rates
oob.error.data <- data.frame(
  Trees=rep(1:nrow(model$err.rate), times=5),
  Type=rep(c("OOB", "C", "C+WBV", "T", "T+WBV"), each=nrow(model$err.rate)),
  Error=c(model$err.rate[,"OOB"],
          model$err.rate[,"C"],
          model$err.rate[,"C+WBV"],
          model$err.rate[,"T"],
          model$err.rate[,"T+WBV"]))

ggplot(data=oob.error.data, aes(x=Trees, y=Error)) +
  geom_line(aes(color=Type)) +
  ggtitle("Plot of error rates Random Forest")

# Check which value for mtry is the best
oob.values <- vector(length=63)
for(i in 1:63) {
  temp.model <- randomForest(as.factor(Group) ~ ., data = DF.groups, ntree=1000, mtry=i)
  oob.values[i] <- temp.model$err.rate[nrow(temp.model$err.rate),1]
}

# Create an MDS-plot to show how the samples are related to each other.
## Converting the proximity matrix into a distance matrix.
distance.matrix <- dist(1-model$proximity)
mds.stuff <- cmdscale(distance.matrix, eig=TRUE, x.ret=TRUE)

## calculate the percentage of variation that each MDS axis accounts for.
mds.var.per <- round(mds.stuff$eig/sum(mds.stuff$eig)*100, 1)

## Make a fancy looking plot that shows the MDS axes and the variation.
mds.values <- mds.stuff$points
mds.data <- data.frame(Sample=rownames(mds.values),
                       X=mds.values[,1],
                       Y=mds.values[,2],
                       Status=DF.groups$Group)

ggplot(data=mds.data, aes(x=X, y=Y, label=Sample)) +
  geom_text(aes(color=Status)) +
  theme_bw() +
  xlab(paste("MDS1 - ", mds.var.per[1], "%", sep="")) +
  ylab(paste("MDS2 - ", mds.var.per[2], "%", sep="")) +
  ggtitle("MDS plot (1-Random Forest Proximities)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values=c("#9ede81", "#588a42", "#7cb1f2", "#1e65bd")) + 
  theme(legend.position = "right")

#########################################################################################################

# PCA Analysis
pca.DF.t3 <- prcomp(DF.t3, scale=TRUE)
summary(pca.DF.t3)
biplot(pca.DF.t3, main="Biplot", scale=0)

# make a scree plot with pecentages of each PC
pca.var <- pca.DF.t3$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
barplot(pca.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")

# Make groups for PCA
Group.pca <- data.frame("Group" = c("C","C","C","C","C","C","T","T","T","T","T","T",
                                    "T+WBV","T+WBV","T+WBV","T+WBV","T+WBV","C+WBV","C+WBV","C+WBV","C+WBV","C+WBV","C+WBV"), 
                        "TB" = c("C","C","C","C","C","C","T","T","T","T","T","T",
                                 "T","T","T","T","T","C","C","C","C","C","C"),
                        "Training" = c("C","C","C","C","C","C","C","C","C","C","C","C",
                                       "WBV","WBV","WBV","WBV","WBV","WBV","WBV","WBV","WBV","WBV","WBV"))

# Combine data to make a plot
DF.groups.p <- cbind(DF.t3, Group.pca)
pca.DF.t3$x
DF.groups.pca <- cbind (DF.groups.p, pca.DF.t3$x[,1:4])

# Plot with ggplot PC1, PC2, PC and PC4 per group
PC1PC2 <- ggplot(DF.groups.pca, aes(PC1, PC2, col = Group, fill = Group)) +
  scale_color_manual(values = c("#9ede81", "#588a42", "#7cb1f2", "#1e65bd"))+
  scale_fill_manual(values = c("#9ede81", "#588a42", "#7cb1f2", "#1e65bd"))+
  stat_ellipse(geom = "polygon", col = "black", alpha = 0.5) +
  geom_point(shape = 21, col = "black") +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  ggtitle("PC1 & PC2")

PC1PC3 <- ggplot(DF.groups.pca, aes(PC1, PC3, col = Group, fill = Group)) +
  scale_color_manual(values = c("#9ede81", "#588a42", "#7cb1f2", "#1e65bd"))+
  scale_fill_manual(values = c("#9ede81", "#588a42", "#7cb1f2", "#1e65bd"))+
  stat_ellipse(geom = "polygon", col = "black", alpha = 0.5) +
  geom_point(shape = 21, col = "black") +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC3 - ", pca.var.per[3], "%", sep="")) +
  ggtitle("PC1 & PC3")

PC1PC4 <- ggplot(DF.groups.pca, aes(PC1, PC4, col = Group, fill = Group)) +
  scale_color_manual(values = c("#9ede81", "#588a42", "#7cb1f2", "#1e65bd"))+
  scale_fill_manual(values = c("#9ede81", "#588a42", "#7cb1f2", "#1e65bd"))+
  stat_ellipse(geom = "polygon", col = "black", alpha = 0.5) +
  geom_point(shape = 21, col = "black") +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC4 - ", pca.var.per[4], "%", sep="")) +
  ggtitle("PC1 & PC4")

PC2PC3 <- ggplot(DF.groups.pca, aes(PC2, PC3, col = Group, fill = Group)) +
  scale_color_manual(values = c("#9ede81", "#588a42", "#7cb1f2", "#1e65bd"))+
  scale_fill_manual(values = c("#9ede81", "#588a42", "#7cb1f2", "#1e65bd"))+
  stat_ellipse(geom = "polygon", col = "black", alpha = 0.5) +
  geom_point(shape = 21, col = "black") +
  xlab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  ylab(paste("PC3 - ", pca.var.per[3], "%", sep="")) +
  ggtitle("PC2 & PC3")

PC2PC4 <- ggplot(DF.groups.pca, aes(PC2, PC4, col = Group, fill = Group)) +
  scale_color_manual(values = c("#9ede81", "#588a42", "#7cb1f2", "#1e65bd"))+
  scale_fill_manual(values = c("#9ede81", "#588a42", "#7cb1f2", "#1e65bd"))+
  stat_ellipse(geom = "polygon", col = "black", alpha = 0.5) +
  geom_point(shape = 21, col = "black") +
  xlab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  ylab(paste("PC4 - ", pca.var.per[4], "%", sep="")) +
  ggtitle("PC2 & PC4")

PC3PC4 <- ggplot(DF.groups.pca, aes(PC3, PC4, col = Group, fill = Group)) +
  scale_color_manual(values = c("#9ede81", "#588a42", "#7cb1f2", "#1e65bd"))+
  scale_fill_manual(values = c("#9ede81", "#588a42", "#7cb1f2", "#1e65bd"))+
  stat_ellipse(geom = "polygon", col = "black", alpha = 0.5) +
  geom_point(shape = 21, col = "black") +
  xlab(paste("PC3 - ", pca.var.per[3], "%", sep="")) +
  ylab(paste("PC4 - ", pca.var.per[4], "%", sep="")) +
  ggtitle("PC3 & PC4")

ggarrange(PC1PC2, PC1PC3, PC1PC4, PC2PC3, PC2PC4, PC3PC4 + rremove("x.text"), 
          labels = c("A", "B", "C", "D", "E", "F"),
          ncol = 2, nrow = 3)
