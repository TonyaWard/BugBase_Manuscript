#### Create Images for Figure 2 ####

#Load HMP tables
HMP <- read.table("BB_HMP/predicted_phenotypes/predictions.txt", 
				sep='\t', 
				head=T, 
				row=1, 
				check=F, 
				comment="")
HMP_Map <- read.table("data/HMP_Map.txt",
				sep='\t',
				head=T,
				row=1,
				check=F,
				comment='')

#Keep only samples in map and OTU, order the samples
intersect_btwn <- intersect(rownames(HMP_Map),rownames(HMP))
hmp_map <- HMP_Map[intersect_btwn,]
hmp <- droplevels(as.data.frame(HMP[intersect_btwn,]))
colnames(hmp) <- colnames(HMP)
hmp_map <- hmp_map[rownames(hmp),]

#Set the order of the groups
groups <- c("Stool,Subgingival_plaque,Supragingival_plaque,Tongue_dorsum")
groups <- strsplit(groups, ",")[[1]]
ix.keep <- hmp_map[,"HMPBODYSUBSITE"] %in% groups
hmp_map <- droplevels(as.data.frame(hmp_map[ix.keep,]))
hmp <- as.data.frame(hmp[ix.keep,])
colnames(hmp) <- colnames(HMP)
hmp_map[,"HMPBODYSUBSITE"] <- factor(hmp_map[,"HMPBODYSUBSITE"], levels = c(groups))

# Load Yats data
GG <- read.table("BB_Yats/predicted_phenotypes/predictions.txt", 
				sep='\t', 
				head=T, 
				row=1, 
				check=F, 
				comment="")
GG_Map <- read.table("data/Yats_Adult_Map.txt",
				sep='\t',
				head=T,
				row=1,
				check=F,
				comment='')
#Keep only samples in map and OTU, order the samples
intersect_btwn <- intersect(rownames(GG_Map),rownames(GG))
gg_map <- GG_Map[intersect_btwn,]
gg <- droplevels(as.data.frame(GG[intersect_btwn,]))
colnames(gg) <- colnames(GG)
gg_map <- gg_map[rownames(gg),]

#Set the order of the groups
groups <- c("Malawi,Venezuela,USA")
groups <- strsplit(groups, ",")[[1]]
ix.keep <- gg_map[,"COUNTRY"] %in% groups
gg_map <- droplevels(as.data.frame(gg_map[ix.keep,]))
gg <- as.data.frame(gg[ix.keep,])
colnames(gg) <- colnames(GG)
gg_map[,"COUNTRY"] <- factor(gg_map[,"COUNTRY"], levels = c(groups))

#Load Yellow Stone data
Yellow_Stone <- read.table("BB_Yellow/predicted_phenotypes/predictions.txt", 
				sep='\t', 
				head=T, 
				row=1, 
				check=F, 
				comment="")
Yellow_Map <- read.table("data/Yellow_Stone_Map.txt",
				sep='\t',
				head=T,
				row=1,
				check=F,
				comment='')
intersect_btwn <- intersect(rownames(Yellow_Map),rownames(Yellow_Stone))
yellow_map <- Yellow_Map[intersect_btwn,]
yellow_stone <- droplevels(as.data.frame(Yellow_Stone[intersect_btwn,]))
colnames(yellow_stone) <- colnames(Yellow_Stone)
yellow_map <- yellow_map[rownames(yellow_stone),]

#Load Soil data
Soil_PH <- read.table("BB_Soil/predicted_phenotypes/predictions.txt", 
				sep='\t', 
				head=T, 
				row=1, 
				check=F, 
				comment="")
Soil_Map <- read.table("data/Soil_Map.txt",
				sep='\t',
				head=T,
				row=1,
				check=F,
				comment='')
intersect_btwn <- intersect(rownames(Soil_Map),rownames(Soil_PH))
soil_map <- Soil_Map[intersect_btwn,]
soil_ph <- droplevels(as.data.frame(Soil_PH[intersect_btwn,]))
colnames(soil_ph) <- colnames(Soil_PH)
soil_map <- soil_map[rownames(soil_ph),]

# Make the linear models.  This is exactly how it is done in BugBase itself
lm.yellow.aer <- lm(yellow_stone[,"Aerobic"] ~ yellow_map[,"temp"])
lm.yellow.fa <- lm(yellow_stone[,"Facultatively_Anaerobic"] ~ yellow_map[,"temp"])
lm.soil.aer <- lm(soil_ph[,"Aerobic"] ~ soil_map[,"ph"])
lm.soil.st <- lm(soil_ph[,"Stress_Tolerant"] ~ soil_map[,"ph"])

# Set colors and load beeswarm library to make plots
library('RColorBrewer')
library('beeswarm')
cols <- sprintf('%s90',brewer.pal(9,'Set1'))
Pal <- colorRampPalette(c(cols[2],cols[1]))
soil_map$Col <- cols[1]
yellow_map$Col <-  cols[2]

# Create smaller versions of the BugBase plots the right dimensions needed in one PDF
# To this file you will add the stats bars, and re-arrange the group labels to fit.
pdf("Figure_2_Images.pdf", width=8, height=7, useDingbats=FALSE)
par(mfrow = c(3, 4))

beeswarm(hmp[,"Anaerobic"] ~ hmp_map[,"HMPBODYSUBSITE"], 
			corral='random', ylim=c(0,1), cex.axis=1.0, tck=-.05, pch=16, col=cols, 
			xlab='', ylab='',cex=1.0, cex.lab=1, las=2)
bxplot(hmp[,"Anaerobic"] ~ hmp_map[,"HMPBODYSUBSITE"], add=TRUE, lwd=0.75)
title(ylab = "Relative Abundance", line = 3, cex.lab = 1)
title(main="Anaerobic", cex.main=1.0, line=0.5)

beeswarm(hmp[,"Forms_Biofilms"] ~ hmp_map[,"HMPBODYSUBSITE"],
			corral='random',ylim=c(0,1), cex.axis=1.0, tck=-.05, pch=16, col=cols, 
			xlab='', ylab='', cex=1.0, cex.lab=1.0, las=2)
bxplot(hmp[,"Forms_Biofilms"] ~ hmp_map[,"HMPBODYSUBSITE"], add=TRUE, lwd=0.75)
title(main="Forms Biofilms", cex.main=1.0, line=0.5)

beeswarm(gg[,"Aerobic"] ~ gg_map[,"COUNTRY"],
			corral='random', cex.axis=1.0, tck=-.05, pch=16, col=cols, xlab='',
			ylab='', cex=1.0, cex.lab=1.0, las=2, ylim=c(0,0.2))
bxplot(gg[,"Aerobic"] ~ gg_map[,"COUNTRY"], add=TRUE, lwd=0.75)
title(ylab = "Relative Abundance", line = 3, cex.lab = 1)
title(main="Aerobic", cex.main=1.0, line=0.5)

beeswarm(gg[,"Facultatively_Anaerobic"] ~ gg_map[,"COUNTRY"],
			corral='random', cex.axis=1.0, tck=-.05, pch=16, col=cols, xlab='',
			ylab='', cex=1.0, cex.lab=1.0, las=2, ylim=c(0,0.2))
bxplot(gg[,"Facultatively_Anaerobic"] ~ gg_map[,"COUNTRY"], add=TRUE, lwd=0.75)
title(main="Facultative Anaerobic", cex.main=1.0, line=0.5)

beeswarm(hmp[,"Contains_Mobile_Elements"] ~ hmp_map[,"HMPBODYSUBSITE"],
			corral='random',ylim=c(0,1), cex.axis=1.0, tck=-.05, pch=16, col=cols, 
			xlab='', ylab='', cex=1.0, cex.lab=1.0, las=2)
bxplot(hmp[,"Contains_Mobile_Elements"] ~ hmp_map[,"HMPBODYSUBSITE"], add=TRUE, lwd=0.75)
title(ylab = "Relative Abundance", line = 3, cex.lab = 1.0)
title(main="Mobile Element Containing", cex.main=1.0, line=0.5)

beeswarm(hmp[,"Pathogenic"] ~ hmp_map[,"HMPBODYSUBSITE"],
			corral='random',ylim=c(0,1), cex.axis=1.0, tck=-.05, pch=16, col=cols, 
			xlab='', ylab='', cex=1.0, cex.lab=1.0, las=2)
bxplot(hmp[,"Pathogenic"] ~ hmp_map[,"HMPBODYSUBSITE"], add=TRUE, lwd=0.75)
title(main="Pathogenic", cex.main=1.0, line=0.5)

beeswarm(gg[,"Gram_Negative"] ~ gg_map[,"COUNTRY"],
			corral='random', cex.axis=1.0,tck=-.05, pch=16,col=cols, xlab='',
			ylab='', cex=1.0, cex.lab=1.0, las=2, ylim=c(0,1))
bxplot(gg[,"Gram_Negative"] ~ gg_map[,"COUNTRY"], add=TRUE, lwd=0.75)
title(ylab = "Relative Abundance", line = 3, cex.lab = 1.0)
title(main="Gram Negative", cex.main=1.0, line=0.5)

beeswarm(gg[,"Pathogenic"] ~ gg_map[,"COUNTRY"],
			corral='random', cex.axis=1.0, tck=-.05, pch=16, col=cols, xlab='', 
			ylab='', cex=1.0, cex.lab=1.0, las=2, ylim=c(0,1))
bxplot(gg[,"Pathogenic"] ~ gg_map[,"COUNTRY"], add=TRUE, lwd=0.75)
title(main="Pathogenic", cex.main=1.0, line=0.5)

plot(yellow_stone[,"Aerobic"] ~ yellow_map[,"temp"],
			col = yellow_map$Col, tck=-.05, pch=16,cex=1.0, xaxt='n', xlab='',
			ylab='Relative Abundance', cex.lab=1.0, cex.axis=1.0, las=1, ylim=c(0,1))
abline(lm.yellow.aer, lty=2, ylim=c(0,1))
axis(1, cex.axis=1.0, padj=-.095)
title(xlab="Temperature", cex.lab=1.0, line=2)
title(main="Aerobic", cex.main=1.0, line=0.5)

plot(yellow_stone[,"Facultatively_Anaerobic"] ~ yellow_map[,"temp"],
			col = yellow_map$Col, tck=-.05,pch=16, cex=1.0, xaxt='n', 
			xlab='', ylab='', cex.lab=1.0, cex.axis=1.0,las=1, ylim=c(0,1))
abline(lm.yellow.fa, lty=2, ylim=c(0,1))
axis(1, cex.axis=1.0, padj=-.095)
title(xlab="Temperature", cex.lab=1.0, line=2)
title(main="Facultative Aerobic", cex.main=1.0, line=0.5)

plot(soil_ph[,"Aerobic"] ~ soil_map[,"ph"],
			col = soil_map$Col, tck=-.05, pch=16, cex=1.0, xaxt='n', xlab='',
			ylab='Relative Abundance', cex.lab=1.0, cex.axis=1.0,las=1, ylim=c(0.2,1))
abline(lm.soil.aer, lty=2, ylim=c(0.2,1))
axis(1, cex.axis=1.0, padj=-.095)
title(xlab="pH", cex.lab=1.0, line=2)
title(main="Aerobic", cex.main=1.0, line=0.5)

plot(soil_ph[,"Stress_Tolerant"] ~ soil_map[,"ph"],
		col = soil_map$Col, tck=-.05,pch=16, cex=1.0, xaxt='n', xlab='', 
		ylab='', cex.lab=1.0, cex.axis=1.0,las=1, ylim=c(0.2,1))
abline(lm.soil.st, lty=2, ylim=c(0.2,1))
axis(1, cex.axis=1.0, padj=-.08)
title(xlab="pH", cex.lab=1.0, line=2)
title(main="Stress_Tolerant", cex.main=1.0, line=0.5)
dev.off()
