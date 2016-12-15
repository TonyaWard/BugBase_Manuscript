
#### Analysis for Figure 3D and Supplemental Figure X ####
#### Compares Bag of Genes (PICRUSt) to BugBase ####

# Run PICRUSt on 16S normalized HMP data used in the Bugbase analysis (BB_HMP/normalized_otus/16s_normalized_otus.txt)

# Load libraries needed
library('beeswarm')
library('RColorBrewer')
library('ggplot2')
library('reshape')
library('gridExtra')
library('grid')
library('lattice')

# Load Tables
# PICRUSt output was generated using the following command:
# predict_metagenomes.py -i 16s_normalized_otus.txt -o predicted_HMP_PICRUSt_pathways.txt -f

picrust <- read.table("data/predicted_HMP_PICRUSt_pathways.txt", 
						sep='\t', 
						skip=1, 
						head=T, 
						row=1, 
						check=F, 
						comment="", 
						quote="")

picrust$KEGG_Description <- NULL # remove the Kegg column that is text
map <- read.table("data/HMP_Map.txt", 
					sep='\t',
					head=T,
					row=1,
					check=F,
					comment='')

# Make the gene abundances relative to all genes present (columns sum to 1)
normalize.picrust <- as.matrix(picrust)
normalize.picrust <- prop.table(normalize.picrust, 2)

# Read in Kegg IDs
kegg_pathways <- read.table("data/kegg_pathways.txt", header=TRUE, sep="\t")

# Subset sample ids into their bodysites
tongue_ids <- rownames(subset(map, HMPBODYSUBSITE == "Tongue_dorsum"))
subging_ids <- rownames(subset(map, HMPBODYSUBSITE == "Subgingival_plaque"))
supging_ids <- rownames(subset(map, HMPBODYSUBSITE == "Supragingival_plaque"))
stool_ids <- rownames(subset(map, HMPBODYSUBSITE == "Stool"))

# Get keggIDs and index them in the PICRUSt table and then make new tables of those IDs per trait
kegg_trypto <- subset(kegg_pathways, Pathway == "Tryptophan_Metabolism")
kegg_benzo <- subset(kegg_pathways, Pathway == "Benzoate_Degradation")
kegg_galact <- subset(kegg_pathways, Pathway == "Galactose_Metabolism")
kegg_gluta <- subset(kegg_pathways, Pathway == "Glutathione_Metabolism")
kegg_starch <- subset(kegg_pathways, Pathway == "Starch_Metabolism")

trypto_ids <- kegg_trypto$Kegg_ID
benzo_ids <- kegg_benzo$Kegg_ID
galact_ids <- kegg_galact$Kegg_ID
gluta_ids <- kegg_gluta$Kegg_ID
starch_ids <- kegg_starch$Kegg_ID

trypto.ix <- rownames(normalize.picrust) %in% trypto_ids
benzo.ix <- rownames(normalize.picrust) %in% benzo_ids
galact.ix <- rownames(normalize.picrust) %in% galact_ids
gluta.ix <- rownames(normalize.picrust) %in% gluta_ids
starch.ix <- rownames(normalize.picrust) %in% starch_ids

picrust_trypto <- droplevels(as.data.frame(normalize.picrust[trypto.ix,]))
picrust_benzo <- droplevels(as.data.frame(normalize.picrust[benzo.ix,]))
picrust_galact <- droplevels(as.data.frame(normalize.picrust[galact.ix,]))
picrust_gluta <- droplevels(as.data.frame(normalize.picrust[gluta.ix,]))
picrust_starch <- droplevels(as.data.frame(normalize.picrust[starch.ix,]))

#Get column sums as new row
picrust_trypto[is.na(picrust_trypto)] <- 0
picrust_trypto["sum",]<- colSums(picrust_trypto, na.rm=TRUE)
picrust_benzo[is.na(picrust_benzo)] <- 0
picrust_benzo["sum",]<- colSums(picrust_benzo, na.rm=TRUE)
picrust_galact[is.na(picrust_galact)] <- 0
picrust_galact["sum",]<- colSums(picrust_galact, na.rm=TRUE)
picrust_gluta[is.na(picrust_gluta)] <- 0
picrust_gluta["sum",]<- colSums(picrust_gluta, na.rm=TRUE)
picrust_starch[is.na(picrust_starch)] <- 0
picrust_starch["sum",]<- colSums(picrust_starch, na.rm=TRUE)

#create master table of all traits
picrust_overall <- picrust_trypto[10,]
rownames(picrust_overall) <- "Tryptophan_Metabolism.txt"
picrust_overall["Benzoate_Degradation.txt",] <- picrust_benzo['sum',]
picrust_overall["Galactose_Metabolism.txt",] <- picrust_galact['sum',]
picrust_overall["Glutathione_Metabolism.txt",] <- picrust_gluta['sum',]
picrust_overall["Starch_Metabolism.txt",] <- picrust_starch['sum',]

picrust_overall <- t(picrust_overall)
picrust_overall[is.na(picrust_overall)] <- 0

intersect_btwn <- intersect(rownames(map),rownames(picrust_overall))

#Keep only samples that intersect between the mapping file and prediction tables
picrust_overall1 <- droplevels(as.data.frame(picrust_overall[intersect_btwn,]))
colnames(picrust_overall1) <- colnames(picrust_overall1)

# read in bug base predictions and the mapping file
traits <- read.table("Fig_3_Pathways/predicted_phenotypes/predictions.txt",
						sep='\t', 
						head=T, 
						row=1, 
						check=F, 
						comment="")
intersect_btwn <- intersect(rownames(map),rownames(traits))
new_map <- map[intersect_btwn,]
new_traits <- droplevels(as.data.frame(traits[intersect_btwn,]))
new_map <- new_map[rownames(new_traits),]

# Subset tables to be body-site specific
stool_picrust <- picrust_overall1[rownames(picrust_overall1) %in% stool_ids,]
tongue_picrust <- picrust_overall1[rownames(picrust_overall1) %in% tongue_ids,]
subg_picrust <- picrust_overall1[rownames(picrust_overall1) %in% subging_ids,]
suprag_picrust <- picrust_overall1[rownames(picrust_overall1) %in% supging_ids,]

stool_bb <- new_traits[rownames(new_traits) %in% stool_ids,]
tongue_bb <- new_traits[rownames(new_traits) %in% tongue_ids,]
subg_bb <- new_traits[rownames(new_traits) %in% subging_ids,]
suprag_bb <- new_traits[rownames(new_traits) %in% supging_ids,]

#Set what to compare
map_column <- "HMPBODYSUBSITE"
#Set order of traits to ensure same comparisons
traits <- colnames(picrust_overall)

print("Making Comparisons, this can take up to 20 minutes")

#Make comparisons for each trait by using 1-250 samples per group
#This loop is for comparing subging plaque to supraging plaque
for(k in 1:length(traits)){
	trait <- traits[k]
	pval_table <- data.frame(matrix(NA, nrow = 6, ncol = 250))
	rownames(pval_table) <- c("bb_mean", "pic_mean", "bb_max", "pic_max", "bb_min", "pic_min")
	colnames(pval_table) <- c(1:250)
	for(i in 1:250){
		picrust_values <- c()
		bb_values <- c()
		for(v in 1:30){ #For each group size, do 30 seperate comparisons using random samples within each group
			sam_sub_bb <- subg_bb[sample(1:nrow(subg_bb), size=i, replace=FALSE),]
			sam_sub_picrust <- subg_picrust[sample(1:nrow(subg_picrust), size=i, replace=FALSE),]
			sam_sup_bb <- suprag_bb[sample(1:nrow(suprag_bb), size=i, replace=FALSE),]
			sam_sup_picrust <- suprag_picrust[sample(1:nrow(suprag_picrust), size=i, replace=FALSE),]
			combine_picr <- rbind(sam_sub_picrust, sam_sup_picrust)
			combine_bbr <- rbind(sam_sub_bb, sam_sup_bb)
			bb_mapr <- new_map[rownames(new_map) %in% rownames(combine_bbr),]
			bb_mapr <- bb_mapr[rownames(combine_bbr),]
			bb_mapr[,map_column] <- factor(bb_mapr[,map_column])
			pic_mapr <- new_map[rownames(new_map) %in% rownames(combine_picr),]
			pic_mapr <- pic_mapr[rownames(combine_picr),]
			pic_mapr[,map_column] <- factor(pic_mapr[,map_column])
			bb_pvalr <- kruskal.test(combine_bbr[,trait] ~ bb_mapr[,map_column])$p.value
			pic_pvalr <- kruskal.test(combine_picr[,trait] ~ pic_mapr[,map_column])$p.value
			picrust_values <- c(picrust_values, pic_pvalr)
			bb_values <- c(bb_values, bb_pvalr)
		}
		pval_table["bb_mean", i] <- mean(bb_values)
		pval_table["pic_mean", i] <- mean(picrust_values)
		pval_table["bb_max", i] <- mean(bb_values) + sd(bb_values)
		pval_table["pic_max", i] <- mean(picrust_values) + sd(picrust_values)
		pval_table["bb_min", i] <- mean(bb_values) - sd(bb_values)
		pval_table["pic_min", i] <- mean(picrust_values) - sd(picrust_values)
		}
		traitx <- paste(trait, "a", sep="") #store the tables with an "a" at the end for 
		assign(traitx, pval_table)	
}

#Make comparisons for each trait by using 1-250 samples per group
#This loop is for comparing stool to tongue
for(k in 1:length(traits)){
	trait <- traits[k]
	pval_table <- data.frame(matrix(NA, nrow = 6, ncol = 250))
	rownames(pval_table) <- c("bb_mean", "pic_mean", "bb_max", "pic_max", "bb_min", "pic_min")
	colnames(pval_table) <- c(1:250)
	for(i in 1:250){
		picrust_values <- c()
		bb_values <- c()
		for(v in 1:30){#For each group size, do 30 seperate comparisons using random samples within each group
			sam_stool_bb <- stool_bb[sample(1:nrow(stool_bb), size=i, replace=FALSE),]
			sam_stool_picrust <- stool_picrust[sample(1:nrow(stool_picrust), size=i, replace=FALSE),]
			sam_tongue_bb <- tongue_bb[sample(1:nrow(tongue_bb), size=i, replace=FALSE),]
			sam_tongue_picrust <- tongue_picrust[sample(1:nrow(tongue_picrust), size=i, replace=FALSE),]
			combine_pic <- rbind(sam_stool_picrust, sam_tongue_picrust)
			combine_bb <- rbind(sam_stool_bb, sam_tongue_bb)
			bb_map <- new_map[rownames(new_map) %in% rownames(combine_bb),]
			bb_map <- bb_map[rownames(combine_bb),]
			bb_map[,map_column] <- factor(bb_map[,map_column])
			pic_map <- new_map[rownames(new_map) %in% rownames(combine_pic),]
			pic_map <- pic_map[rownames(combine_pic),]
			pic_map[,map_column] <- factor(pic_map[,map_column])
			bb_pval <- kruskal.test(combine_bb[,trait] ~ bb_map[,map_column])$p.value
			pic_pval <- kruskal.test(combine_pic[,trait] ~ pic_map[,map_column])$p.value
			picrust_values <- c(picrust_values, pic_pval)
			bb_values <- c(bb_values, bb_pval)
		}
		pval_table["bb_mean", i] <- mean(bb_values)
		pval_table["pic_mean", i] <- mean(picrust_values)
		pval_table["bb_max", i] <- mean(bb_values) + sd(bb_values)
		pval_table["pic_max", i] <- mean(picrust_values) + sd(picrust_values)
		pval_table["bb_min", i] <- mean(bb_values) - sd(bb_values)
		pval_table["pic_min", i] <- mean(picrust_values) - sd(picrust_values)
		}
		assign(trait, pval_table) #these tables don't have the 'a' at the end
}

#Transpose
Benzoate_Degradation.txt <- t(Benzoate_Degradation.txt) #stool/tongue
Benzoate_Degradation.txta <- t(Benzoate_Degradation.txta) #plaque/plaque
Galactose_Metabolism.txt <- t(Galactose_Metabolism.txt)
Galactose_Metabolism.txta <- t(Galactose_Metabolism.txta)
Glutathione_Metabolism.txt <- t(Glutathione_Metabolism.txt)
Glutathione_Metabolism.txta <- t(Glutathione_Metabolism.txta)
Tryptophan_Metabolism.txt <- t(Tryptophan_Metabolism.txt)
Tryptophan_Metabolism.txta <- t(Tryptophan_Metabolism.txta)
Starch_Metabolism.txt <- t(Starch_Metabolism.txt)
Starch_Metabolism.txta <- t(Starch_Metabolism.txta)

#Make a data frame
Benzoate_Degradation.txt <- as.data.frame(Benzoate_Degradation.txt)
Benzoate_Degradation.txta <- as.data.frame(Benzoate_Degradation.txta)
Galactose_Metabolism.txt <- as.data.frame(Galactose_Metabolism.txt)
Galactose_Metabolism.txta <- as.data.frame(Galactose_Metabolism.txta)
Glutathione_Metabolism.txt <- as.data.frame(Glutathione_Metabolism.txt)
Glutathione_Metabolism.txta <- as.data.frame(Glutathione_Metabolism.txta)
Tryptophan_Metabolism.txt <- as.data.frame(Tryptophan_Metabolism.txt)
Tryptophan_Metabolism.txta <- as.data.frame(Tryptophan_Metabolism.txta)
Starch_Metabolism.txt <- as.data.frame(Starch_Metabolism.txt)
Starch_Metabolism.txta <- as.data.frame(Starch_Metabolism.txta)

#Add n as a column
Benzoate_Degradation.txt$subsample <- seq.int(nrow(Benzoate_Degradation.txt))
Benzoate_Degradation.txta$subsample <- seq.int(nrow(Benzoate_Degradation.txta))
Galactose_Metabolism.txt$subsample <- seq.int(nrow(Galactose_Metabolism.txt))
Galactose_Metabolism.txta$subsample <- seq.int(nrow(Galactose_Metabolism.txta))
Glutathione_Metabolism.txt$subsample <- seq.int(nrow(Glutathione_Metabolism.txt))
Glutathione_Metabolism.txta$subsample <- seq.int(nrow(Glutathione_Metabolism.txta))
Tryptophan_Metabolism.txt$subsample <- seq.int(nrow(Tryptophan_Metabolism.txt))
Tryptophan_Metabolism.txta$subsample <- seq.int(nrow(Tryptophan_Metabolism.txta))
Starch_Metabolism.txt$subsample <- seq.int(nrow(Starch_Metabolism.txt))
Starch_Metabolism.txta$subsample <- seq.int(nrow(Starch_Metabolism.txta))

#Rename the traits
Benzoate_Degradation.txt$pathway <- "Benzoate_Degradation"
Benzoate_Degradation.txta$pathway <- "Benzoate_Degradation"
Galactose_Metabolism.txt$pathway <- "Galactose_Metabolism"
Galactose_Metabolism.txta$pathway <- "Galactose_Metabolism"
Glutathione_Metabolism.txt$pathway <- "Glutathione_Metabolism"
Glutathione_Metabolism.txta$pathway <- "Glutathione_Metabolism"
tophan_Metabolism.txt$pathway <- "Tryptophan_Metabolism"
Tryptophan_Metabolism.txta$pathway <- "Tryptophan_Metabolism"
Starch_Metabolism.txt$pathway <- "Starch_Metabolism"
Starch_Metabolism.txta$pathway <- "Starch_Metabolism"

#Smooth the data
Tryptophan_Metabolism.txt[is.na(Tryptophan_Metabolism.txt)] <- 0
mean_lineb <- predict(loess(Tryptophan_Metabolism.txt$bb_mean ~ Tryptophan_Metabolism.txt$subsample))
top_lineb <- predict(loess(Tryptophan_Metabolism.txt$bb_max ~ Tryptophan_Metabolism.txt$subsample))
bottom_lineb <- predict(loess(Tryptophan_Metabolism.txt$bb_min ~ 
Tryptophan_Metabolism.txt$subsample))
mean_linep <- predict(loess(Tryptophan_Metabolism.txt$pic_mean ~ Tryptophan_Metabolism.txt$subsample))
top_linep <- predict(loess(Tryptophan_Metabolism.txt$pic_max ~ Tryptophan_Metabolism.txt$subsample))
bottom_linep <- predict(loess(Tryptophan_Metabolism.txt$pic_min ~ Tryptophan_Metabolism.txt$subsample))
tryp_df <- data.frame(Tryptophan_Metabolism.txt$subsample,top_lineb,bottom_lineb,mean_lineb,top_linep,bottom_linep,mean_linep)

Galactose_Metabolism.txt[is.na(Galactose_Metabolism.txt)] <- 0
mean_lineb <- predict(loess(Galactose_Metabolism.txt$bb_mean ~ Galactose_Metabolism.txt$subsample))
top_lineb <- predict(loess(Galactose_Metabolism.txt$bb_max ~ Galactose_Metabolism.txt$subsample))
bottom_lineb <- predict(loess(Galactose_Metabolism.txt$bb_min ~ 
Galactose_Metabolism.txt$subsample))
mean_linep <- predict(loess(Galactose_Metabolism.txt$pic_mean ~ Galactose_Metabolism.txt$subsample))
top_linep <- predict(loess(Galactose_Metabolism.txt$pic_max ~ Galactose_Metabolism.txt$subsample))
bottom_linep <- predict(loess(Galactose_Metabolism.txt$pic_min ~ Galactose_Metabolism.txt$subsample))
gal_df <- data.frame(Galactose_Metabolism.txt$subsample,top_lineb,bottom_lineb,mean_lineb,top_linep,bottom_linep,mean_linep)

Benzoate_Degradation.txt[is.na(Benzoate_Degradation.txt)] <- 0
mean_lineb <- predict(loess(Benzoate_Degradation.txt$bb_mean ~ Benzoate_Degradation.txt$subsample))
top_lineb <- predict(loess(Benzoate_Degradation.txt$bb_max ~ Benzoate_Degradation.txt$subsample))
bottom_lineb <- predict(loess(Benzoate_Degradation.txt$bb_min ~ 
Benzoate_Degradation.txt$subsample))
mean_linep <- predict(loess(Benzoate_Degradation.txt$pic_mean ~ Benzoate_Degradation.txt$subsample))
top_linep <- predict(loess(Benzoate_Degradation.txt$pic_max ~ Benzoate_Degradation.txt$subsample))
bottom_linep <- predict(loess(Benzoate_Degradation.txt$pic_min ~ Benzoate_Degradation.txt$subsample))
ben_df <- data.frame(Benzoate_Degradation.txt$subsample,top_lineb,bottom_lineb,mean_lineb,top_linep,bottom_linep,mean_linep)

Starch_Metabolism.txt[is.na(Starch_Metabolism.txt)] <- 0
mean_lineb <- predict(loess(Starch_Metabolism.txt$bb_mean ~ Starch_Metabolism.txt$subsample))
top_lineb <- predict(loess(Starch_Metabolism.txt$bb_max ~ Starch_Metabolism.txt$subsample))
bottom_lineb <- predict(loess(Starch_Metabolism.txt$bb_min ~ 
Starch_Metabolism.txt$subsample))
mean_linep <- predict(loess(Starch_Metabolism.txt$pic_mean ~ Starch_Metabolism.txt$subsample))
top_linep <- predict(loess(Starch_Metabolism.txt$pic_max ~ Starch_Metabolism.txt$subsample))
bottom_linep <- predict(loess(Starch_Metabolism.txt$pic_min ~ Starch_Metabolism.txt$subsample))
starch_df <- data.frame(Starch_Metabolism.txt$subsample,top_lineb,bottom_lineb,mean_lineb,top_linep,bottom_linep,mean_linep)

Glutathione_Metabolism.txt[is.na(Glutathione_Metabolism.txt)] <- 0
mean_lineb <- predict(loess(Glutathione_Metabolism.txt$bb_mean ~ Glutathione_Metabolism.txt$subsample))
top_lineb <- predict(loess(Glutathione_Metabolism.txt$bb_max ~ Glutathione_Metabolism.txt$subsample))
bottom_lineb <- predict(loess(Glutathione_Metabolism.txt$bb_min ~ 
Glutathione_Metabolism.txt$subsample))
mean_linep <- predict(loess(Glutathione_Metabolism.txt$pic_mean ~ Glutathione_Metabolism.txt$subsample))
top_linep <- predict(loess(Glutathione_Metabolism.txt$pic_max ~ Glutathione_Metabolism.txt$subsample))
bottom_linep <- predict(loess(Glutathione_Metabolism.txt$pic_min ~ Glutathione_Metabolism.txt$subsample))
gluta_df <- data.frame(Glutathione_Metabolism.txt$subsample,top_lineb,bottom_lineb,mean_lineb,top_linep,bottom_linep,mean_linep)

Tryptophan_Metabolism.txta[is.na(Tryptophan_Metabolism.txta)] <- 0
mean_lineb <- predict(loess(Tryptophan_Metabolism.txta$bb_mean ~ Tryptophan_Metabolism.txta$subsample))
top_lineb <- predict(loess(Tryptophan_Metabolism.txta$bb_max ~ Tryptophan_Metabolism.txta$subsample))
bottom_lineb <- predict(loess(Tryptophan_Metabolism.txta$bb_min ~ 
Tryptophan_Metabolism.txta$subsample))
mean_linep <- predict(loess(Tryptophan_Metabolism.txta$pic_mean ~ Tryptophan_Metabolism.txta$subsample))
top_linep <- predict(loess(Tryptophan_Metabolism.txta$pic_max ~ Tryptophan_Metabolism.txta$subsample))
bottom_linep <- predict(loess(Tryptophan_Metabolism.txta$pic_min ~ Tryptophan_Metabolism.txta$subsample))
tryp_dfa <- data.frame(Tryptophan_Metabolism.txta$subsample,top_lineb,bottom_lineb,mean_lineb,top_linep,bottom_linep,mean_linep)

Galactose_Metabolism.txta[is.na(Galactose_Metabolism.txta)] <- 0
mean_lineb <- predict(loess(Galactose_Metabolism.txta$bb_mean ~ Galactose_Metabolism.txta$subsample))
top_lineb <- predict(loess(Galactose_Metabolism.txta$bb_max ~ Galactose_Metabolism.txta$subsample))
bottom_lineb <- predict(loess(Galactose_Metabolism.txta$bb_min ~ 
Galactose_Metabolism.txta$subsample))
mean_linep <- predict(loess(Galactose_Metabolism.txta$pic_mean ~ Galactose_Metabolism.txta$subsample))
top_linep <- predict(loess(Galactose_Metabolism.txta$pic_max ~ Galactose_Metabolism.txta$subsample))
bottom_linep <- predict(loess(Galactose_Metabolism.txta$pic_min ~ Galactose_Metabolism.txta$subsample))
gal_dfa <- data.frame(Galactose_Metabolism.txta$subsample,top_lineb,bottom_lineb,mean_lineb,top_linep,bottom_linep,mean_linep)

Benzoate_Degradation.txta[is.na(Benzoate_Degradation.txta)] <- 0
mean_lineb <- predict(loess(Benzoate_Degradation.txta$bb_mean ~ Benzoate_Degradation.txta$subsample))
top_lineb <- predict(loess(Benzoate_Degradation.txta$bb_max ~ Benzoate_Degradation.txta$subsample))
bottom_lineb <- predict(loess(Benzoate_Degradation.txta$bb_min ~ 
Benzoate_Degradation.txta$subsample))
mean_linep <- predict(loess(Benzoate_Degradation.txta$pic_mean ~ Benzoate_Degradation.txta$subsample))
top_linep <- predict(loess(Benzoate_Degradation.txta$pic_max ~ Benzoate_Degradation.txta$subsample))
bottom_linep <- predict(loess(Benzoate_Degradation.txta$pic_min ~ Benzoate_Degradation.txta$subsample))
ben_dfa <- data.frame(Benzoate_Degradation.txta$subsample,top_lineb,bottom_lineb,mean_lineb,top_linep,bottom_linep,mean_linep)

Starch_Metabolism.txta[is.na(Starch_Metabolism.txta)] <- 0
mean_lineb <- predict(loess(Starch_Metabolism.txta$bb_mean ~ Starch_Metabolism.txta$subsample))
top_lineb <- predict(loess(Starch_Metabolism.txta$bb_max ~ Starch_Metabolism.txta$subsample))
bottom_lineb <- predict(loess(Starch_Metabolism.txta$bb_min ~ 
Starch_Metabolism.txta$subsample))
mean_linep <- predict(loess(Starch_Metabolism.txta$pic_mean ~ Starch_Metabolism.txta$subsample))
top_linep <- predict(loess(Starch_Metabolism.txta$pic_max ~ Starch_Metabolism.txta$subsample))
bottom_linep <- predict(loess(Starch_Metabolism.txta$pic_min ~ Starch_Metabolism.txta$subsample))
starch_dfa <- data.frame(Starch_Metabolism.txta$subsample,top_lineb,bottom_lineb,mean_lineb,top_linep,bottom_linep,mean_linep)

Glutathione_Metabolism.txta[is.na(Glutathione_Metabolism.txta)] <- 0
mean_lineb <- predict(loess(Glutathione_Metabolism.txta$bb_mean ~ Glutathione_Metabolism.txta$subsample))
top_lineb <- predict(loess(Glutathione_Metabolism.txta$bb_max ~ Glutathione_Metabolism.txta$subsample))
bottom_lineb <- predict(loess(Glutathione_Metabolism.txta$bb_min ~ 
Glutathione_Metabolism.txta$subsample))
mean_linep <- predict(loess(Glutathione_Metabolism.txta$pic_mean ~ Glutathione_Metabolism.txta$subsample))
top_linep <- predict(loess(Glutathione_Metabolism.txta$pic_max ~ Glutathione_Metabolism.txta$subsample))
bottom_linep <- predict(loess(Glutathione_Metabolism.txta$pic_min ~ Glutathione_Metabolism.txta$subsample))
gluta_dfa <- data.frame(Glutathione_Metabolism.txta$subsample,top_lineb,bottom_lineb,mean_lineb,top_linep,bottom_linep,mean_linep)

#Set colors
cols <- brewer.pal(9,'Set1')

#gg-plot it
tryp_plot <- ggplot(tryp_df, aes(x = Tryptophan_Metabolism.txt.subsample)) + 
				theme_bw() + 
				scale_color_brewer(palette="Set1") + 
				labs(x="",y="p-value") + 
				theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), panel.border=element_blank()) + 
				theme(legend.key= element_blank(), legend.title = element_blank()) + 
				guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5, override.aes = list(colour = NULL))) + 
				ggtitle("Tryptophan Metabolism") + 
				geom_line(aes(x=Tryptophan_Metabolism.txt.subsample, y=mean_lineb),color=cols[2], size=1.3)+ 
				geom_ribbon(aes(ymin=bottom_lineb, ymax=top_lineb, linetype=NA), fill=cols[2], alpha=0.3) + 
				geom_line(aes(x=Tryptophan_Metabolism.txt.subsample, y=mean_linep), color=cols[1], size=1.3)+ 
				geom_ribbon(aes(ymin=bottom_linep, ymax=top_linep, linetype=NA), fill=cols[1], alpha=0.3) 

tryp_plota <- ggplot(tryp_dfa, aes(x = Tryptophan_Metabolism.txta.subsample)) + 
				theme_bw() + scale_color_brewer(palette="Set1") + 
				labs(x="",y="p-value") + 
				theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), panel.border=element_blank()) + 
				theme(legend.key= element_blank(), legend.title = element_blank()) + 
				guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5, override.aes = list(colour = NULL))) + 
				ggtitle("Tryptophan Metabolism") + 
				geom_line(aes(x=Tryptophan_Metabolism.txta.subsample, y=mean_lineb),color=cols[2], size=1.3)+ 
				geom_ribbon(aes(ymin=bottom_lineb, ymax=top_lineb, linetype=NA), fill=cols[2], alpha=0.3) + 
				geom_line(aes(x=Tryptophan_Metabolism.txta.subsample, y=mean_linep), color=cols[1], size=1.3)+ 
				geom_ribbon(aes(ymin=bottom_linep, ymax=top_linep, linetype=NA), fill=cols[1], alpha=0.3) 

gal_plot <- ggplot(gal_df, aes(x = Galactose_Metabolism.txt.subsample)) + 
				theme_bw() + 
				scale_color_brewer(palette="Set1") + 
        labs(x="",y="p-value") + 
				theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), panel.border=element_blank()) + 
				theme(legend.key= element_blank(), legend.title = element_blank()) + 
				guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5, override.aes = list(colour = NULL))) + 
				ggtitle("Galactose Metabolism") + 
				geom_line(aes(x=Galactose_Metabolism.txt.subsample, y=mean_lineb),color=cols[2], size=1.3)+ 
				geom_ribbon(aes(ymin=bottom_lineb, ymax=top_lineb, linetype=NA), fill=cols[2], alpha=0.3) + 
				geom_line(aes(x=Galactose_Metabolism.txt.subsample, y=mean_linep), color=cols[1], size=1.3)+ 
				geom_ribbon(aes(ymin=bottom_linep, ymax=top_linep, linetype=NA), fill=cols[1], alpha=0.3) 

gal_plota <- ggplot(gal_dfa, aes(x = Galactose_Metabolism.txta.subsample)) + 
				theme_bw() + 
				scale_color_brewer(palette="Set1") + labs(x="",y="p-value") + 
				theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), panel.border=element_blank()) + 
				theme(legend.key= element_blank(), legend.title = element_blank()) + 
				guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5, override.aes = list(colour = NULL)))+
				ggtitle("Galactose Metabolism") + 
				geom_line(aes(x=Galactose_Metabolism.txta.subsample, y=mean_lineb),color=cols[2], size=1.3)+ 
				geom_ribbon(aes(ymin=bottom_lineb, ymax=top_lineb, linetype=NA), fill=cols[2], alpha=0.3) + 
				geom_line(aes(x=Galactose_Metabolism.txta.subsample, y=mean_linep), color=cols[1], size=1.3)+ 
				geom_ribbon(aes(ymin=bottom_linep, ymax=top_linep, linetype=NA), fill=cols[1], alpha=0.3) 

ben_plot <- ggplot(ben_df, aes(x = Benzoate_Degradation.txt.subsample)) + 
				theme_bw() + 
				scale_color_brewer(palette="Set1") + labs(x="",y="p-value") + 
				theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), panel.border=element_blank()) + 
				theme(legend.key= element_blank(), legend.title = element_blank()) + 
				guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5, override.aes = list(colour = NULL))) + 
				ggtitle("Benzoate Degradation") + 
				geom_line(aes(x=Benzoate_Degradation.txt.subsample, y=mean_lineb),color=cols[2], size=1.3)+ 
				geom_ribbon(aes(ymin=bottom_lineb, ymax=top_lineb, linetype=NA), fill=cols[2], alpha=0.3) + 
				geom_line(aes(x=Benzoate_Degradation.txt.subsample, y=mean_linep), color=cols[1], size=1.3)+ 
				geom_ribbon(aes(ymin=bottom_linep, ymax=top_linep, linetype=NA), fill=cols[1], alpha=0.3) 

ben_plota <- ggplot(ben_dfa, aes(x = Benzoate_Degradation.txta.subsample)) + 
				theme_bw() + 
				scale_color_brewer(palette="Set1") + labs(x="",y="p-value") + 
				theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), panel.border=element_blank()) + 
				theme(legend.key= element_blank(), legend.title = element_blank()) + 
				guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5, override.aes = list(colour = NULL))) + 
				ggtitle("Benzoate Degradation") + 
				geom_line(aes(x=Benzoate_Degradation.txta.subsample, y=mean_lineb),color=cols[2], size=1.3)+ 
				geom_ribbon(aes(ymin=bottom_lineb, ymax=top_lineb, linetype=NA), fill=cols[2], alpha=0.3) + 
				geom_line(aes(x=Benzoate_Degradation.txta.subsample, y=mean_linep), color=cols[1], size=1.3)+ 
				geom_ribbon(aes(ymin=bottom_linep, ymax=top_linep, linetype=NA), fill=cols[1], alpha=0.3) 

starch_plot <- ggplot(starch_df, aes(x = Starch_Metabolism.txt.subsample)) + 
				theme_bw() + 
				scale_color_brewer(palette="Set1") + 
				labs(x="",y="p-value") + 
				theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), panel.border=element_blank()) + 
				theme(legend.key= element_blank(), legend.title = element_blank()) + 
				guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5, override.aes = list(colour = NULL))) + 
				ggtitle("Starch Metabolism") + 
				geom_line(aes(x=Starch_Metabolism.txt.subsample, y=mean_lineb),color=cols[2], size=1.3)+ 
				geom_ribbon(aes(ymin=bottom_lineb, ymax=top_lineb, linetype=NA), fill=cols[2], alpha=0.3) + 
				geom_line(aes(x=Starch_Metabolism.txt.subsample, y=mean_linep), color=cols[1], size=1.3)+ 
				geom_ribbon(aes(ymin=bottom_linep, ymax=top_linep, linetype=NA), fill=cols[1], alpha=0.3) 

starch_plota <- ggplot(starch_dfa, aes(x = Starch_Metabolism.txta.subsample)) + 
				theme_bw() + 
				scale_color_brewer(palette="Set1") + 
				labs(x="",y="p-value") + 
				theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), panel.border=element_blank()) + 
				theme(legend.key= element_blank(), legend.title = element_blank()) + 
				guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5, override.aes = list(colour = NULL))) + 
				ggtitle("Starch Metabolism") + 
				geom_line(aes(x=Starch_Metabolism.txta.subsample, y=mean_lineb),color=cols[2], size=1.3)+ 
				geom_ribbon(aes(ymin=bottom_lineb, ymax=top_lineb, linetype=NA), fill=cols[2], alpha=0.3) + 
				geom_line(aes(x=Starch_Metabolism.txta.subsample, y=mean_linep), color=cols[1], size=1.3)+ 
				geom_ribbon(aes(ymin=bottom_linep, ymax=top_linep, linetype=NA), fill=cols[1], alpha=0.3)

gluta_plot <- ggplot(gluta_df, aes(x = Glutathione_Metabolism.txt.subsample)) + 
				theme_bw() + 
				scale_color_brewer(palette="Set1") + 
				labs(x="",y="p-value") + 
				theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), panel.border=element_blank()) + 
				theme(legend.key= element_blank(), legend.title = element_blank()) + 
				guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5, override.aes = list(colour = NULL))) + 
				ggtitle("Glutathione Metabolism") + 
				geom_line(aes(x=Glutathione_Metabolism.txt.subsample, y=mean_lineb),color=cols[2], size=1.3)+ 
				geom_ribbon(aes(ymin=bottom_lineb, ymax=top_lineb, linetype=NA), fill=cols[2], alpha=0.3) + 
				geom_line(aes(x=Glutathione_Metabolism.txt.subsample, y=mean_linep), color=cols[1], size=1.3)+ 
				geom_ribbon(aes(ymin=bottom_linep, ymax=top_linep, linetype=NA), fill=cols[1], alpha=0.3) 

gluta_plota <- ggplot(gluta_dfa, aes(x = Glutathione_Metabolism.txta.subsample)) + 
				theme_bw() + 
				scale_color_brewer(palette="Set1") + 
				labs(x="",y="p-value") + 
				theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), panel.border=element_blank()) + 
				theme(legend.key= element_blank(), legend.title = element_blank()) + 
				guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5, override.aes = list(colour = NULL))) + 
				ggtitle("Glutathione Metabolism") + 
				geom_line(aes(x=Glutathione_Metabolism.txta.subsample, y=mean_lineb),color=cols[2], size=1.3)+ 
				geom_ribbon(aes(ymin=bottom_lineb, ymax=top_lineb, linetype=NA), fill=cols[2], alpha=0.3) + 
				geom_line(aes(x=Glutathione_Metabolism.txta.subsample, y=mean_linep), color=cols[1], size=1.3)+ 
				geom_ribbon(aes(ymin=bottom_linep, ymax=top_linep, linetype=NA), fill=cols[1], alpha=0.3)

#Make plot all together, save as pdf
pdf("Figure_3_Pathways.pdf", width = 8, height = 12, useDingbats=FALSE)
grid.arrange(tryp_plot, tryp_plota, gal_plot, gal_plota, ben_plot,ben_plota, starch_plot, starch_plota, gluta_plot, gluta_plota, ncol=2)
dev.off()

