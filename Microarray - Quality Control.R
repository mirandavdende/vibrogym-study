# Quality control Vibrogym Experiment

setwd("your directory")

# QC checks
# set colour palette
cols <- brewer.pal(8, "Set1")
# plot a boxplot of unnormalised intensity values
jpeg('boxplot1.jpg')
boxplot(affy.data, col=cols);dev.off()
# plot a boxplot of normalised intensity values, affyPLM is required to interrogate celfiles.gcrma
jpeg('boxplot2.jpg')
boxplot(x.norm, col=cols);dev.off()
# the boxplots are somewhat skewed by the normalisation algorithm
# and it is often more informative to look at density plots
# Plot a density vs log intensity histogram for the unnormalised data
jpeg('hist1.jpg')
hist(affy.data, col=cols);dev.off()
# Plot a density vs log intensity histogram for the normalised data
jpeg('hist2.jpg')
hist(pset2eset(x.norm), col=cols);dev.off()

# Create image of arrays:
# Weights
jpeg('weightsE.jpg',width = 1200, height = 1200, units = "px", pointsize = 12,quality = 100)
op = par(mfrow = c(4,6))
for (i in 1:23){image(x.norm,which=i)};dev.off()
jpeg('weightsS.jpg',width = 1200, height = 1200, units = "px", pointsize = 12,quality = 100)
op = par(mfrow = c(4,6))
for (i in 24:46){image(x.norm,which=i)};dev.off()
jpeg('weightsH.jpg',width = 1200, height = 1200, units = "px", pointsize = 12,quality = 100)
op = par(mfrow = c(4,6))
for (i in 47:58){image(x.norm,which=i)};dev.off()

# Residuals
jpeg('residualsE.jpg',width = 1200, height = 1200, units = "px", pointsize = 12,quality = 100)
op = par(mfrow = c(4,6))
for (i in 1:23){image(x.norm,type='resids',which=i)};dev.off()
jpeg('residualsS.jpg',width = 1200, height = 1200, units = "px", pointsize = 12,quality = 100)
op = par(mfrow = c(4,6))
for (i in 24:46){image(x.norm,type='resids',which=i)};dev.off()
jpeg('residualsH.jpg',width = 1200, height = 1200, units = "px", pointsize = 12,quality = 100)
op = par(mfrow = c(4,6))
for (i in 47:58){image(x.norm,type='resids',which=i)};dev.off()

# Significance
jpeg('significanceE.jpg',width = 1200, height = 1200, units = "px", pointsize = 12,quality = 100)
op = par(mfrow = c(4,6))
for (i in 1:23){image(x.norm,type='sign.resids',which=i)};dev.off()
jpeg('significanceS.jpg',width = 1200, height = 1200, units = "px", pointsize = 12,quality = 100)
op = par(mfrow = c(4,6))
for (i in 24:46){image(x.norm,type='sign.resids',which=i)};dev.off()
jpeg('significanceH.jpg',width = 1200, height = 1200, units = "px", pointsize = 12,quality = 100)
op = par(mfrow = c(4,6))
for (i in 47:58){image(x.norm,type='sign.resids',which=i)};dev.off()

# PCA plot
jpeg('pca.jpg',width = 1200, height = 1200, units = "px", pointsize = 12,quality = 100)
op = par(mfrow= c(1,1))
data.matrix=exprs(pset2eset(x.norm))
color=c('red','orange','cyan','blue','cyan','red',
        'blue','cyan','blue','orange','red','blue',
        'red','cyan','orange','orange','cyan','blue',
        'red','orange','cyan','red','blue',
        'yellow','pink','purple','green','purple','yellow',
        'green','purple','green','pink','yellow','green',
        'yellow','purple','pink','pink','purple','green',
        'yellow','pink','purple','yellow','green',
        'black','brown','gold','grey','gold','black',
        'brown','black','brown','grey','gold','grey')
data.PC = prcomp(t(data.matrix),scale.=TRUE)
plot(data.PC$x,col=color)
text(data.PC$x,labels=samples, pos= 3);
dev.off()
# affyPLM also provides more informative boxplots
# RLE (Relative Log Expression) plots should have
# values close to zero.
jpeg('RLE.jpg')
RLE(x.norm, main="RLE");dev.off()
# We can also use NUSE (Normalised Unscaled Standard Errors).
# The median standard error should be 1 for most genes.
jpeg('NUSE.jpg')
NUSE(x.norm, main="NUSE");dev.off()

