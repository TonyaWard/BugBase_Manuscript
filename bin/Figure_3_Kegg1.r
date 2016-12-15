# Comparing PICRUSt (Bag of Genes) to BugBase
# Run PICRUSt on HMP data using the 16S normalized HMP data used in the Bugbase analysis

library('beeswarm')
library('RColorBrewer')
library('ggplot2')
library('reshape')
library('gridExtra')
library('grid')
library('lattice')

#Read in predictions and maps
picrust <- read.table("data/HMP_PICRUSt.txt", 
                        sep='\t',
                        skip=1,
                        head=T,
                        row=1,
                        check=F,
                        comment="",
                        quote="")
picrust <- picrust[,1:(ncol(picrust)-1)]
picrust$KO_ID <- rownames(picrust)

KO_map <- read.table("data/ko-module-annotations.txt", 
                      sep='\t', 
                      check=F, 
                      comment="", 
                      quote="", 
                      as.is = T)
KO_map$V5 <- gsub('\\"', "", KO_map$V5)

new_table <- matrix(ncol=1512)
colnames(new_table) <- colnames(picrust)

print("Aggregating to KEGG modules, this can take up to 30 minutes")

for(i in 1:nrow(KO_map)){
  KO_ID <- KO_map[i,"V1"]
  MO_ID <- KO_map[i,"V5"]
  for(k in 1:nrow(picrust)){
    if(KO_ID == picrust[k,"KO_ID"]){
      new_table <- rbind(new_table,picrust[k,])
      new_table[nrow(new_table),"KO_ID"] <- MO_ID
    }
  }
}

picrust_m3 <- aggregate(new_table[,1:(ncol(new_table)-1)], by=list(new_table$KO_ID), FUN=sum)
colnames(picrust_m3)[1] <- "module_id"

write.table(picrust_m3, "data/picrust_modules3.txt", 
              quote=FALSE, 
              sep="\t", 
              col.names=NA) 

#picrust_m3 <- read.table("picrust_modules3.txt", sep="\t", head=T, row=1, check=F, comment="", quote="")

BugBase <- read.table("BB_HMP_Kegg/predicted_phenotypes/predictions.txt", 
                        sep="\t", 
                        head=T, 
                        row=1, 
                        check=F, 
                        comment="", 
                        quote="")
BB_map <- read.table("data/Kegg_ID_Map.txt", 
                        sep="\t",
                        head=T,
                        check=F, 
                        comment="", 
                        quote="")
BB_map$KEGG_Label <- as.character(BB_map$KEGG_Label)
colnames(BugBase) <- BB_map$KEGG_Label

#Reorder Picrust to match KO_map
rownames(picrust_m3) <- picrust_m3$module_id
picrust_m3 <- picrust_m3[,2:ncol(picrust_m3)]
picrust_m3 <- as.data.frame(t(picrust_m3))

#Convert to relative abundance to match BB
picrust_m3 <- sweep(picrust_m3,2,colSums(picrust_m3),`/`)
picrust_m3[is.na(picrust_m3)] <- 0

order_Ids <- intersect(colnames(picrust_m3),colnames(BugBase))
picrust_m3 <- picrust_m3[,order_Ids]
BugBase <- BugBase[,order_Ids]

#Load in sample map and subset tables for each site
map <- read.table("data/HMP_Map.txt",
                    sep='\t',head=T,
                    row=1,
                    check=F,
                    comment='')
order_samples <- rownames(map)
picrust_m3 <- picrust_m3[order_samples,]
BugBase <- BugBase[order_samples,]

tongue_ids <- rownames(subset(map, HMPBODYSUBSITE == "Tongue_dorsum"))
subging_ids <- rownames(subset(map, HMPBODYSUBSITE == "Subgingival_plaque"))
supging_ids <- rownames(subset(map, HMPBODYSUBSITE == "Supragingival_plaque"))
stool_ids <- rownames(subset(map, HMPBODYSUBSITE == "Stool"))

stool_picrust <- picrust_m3[rownames(picrust_m3) %in% stool_ids,]
tongue_picrust <- picrust_m3[rownames(picrust_m3) %in% tongue_ids,]
subg_picrust <- picrust_m3[rownames(picrust_m3) %in% subging_ids,]
suprag_picrust <- picrust_m3[rownames(picrust_m3) %in% supging_ids,]

stool_bb <- BugBase[rownames(BugBase) %in% stool_ids,]
tongue_bb <- BugBase[rownames(BugBase) %in% tongue_ids,]
subg_bb <- BugBase[rownames(BugBase) %in% subging_ids,]
suprag_bb <- BugBase[rownames(BugBase) %in% supging_ids,]

map_column <- "HMPBODYSUBSITE"
traits <- colnames(BugBase)

##Pval tables: each row is the mean BB pval and mean Picrust pval for the 30n subsample at that depth
##table wil be 250 (columns) by 470 (rows)
##Plot each mean as a line with transparency
##

#One table for BB pvals and one for picrust
#Compare sub and supra plaque first
pval_table <- data.frame(matrix(NA, nrow = 470, ncol = 5))
sizes <- c(5, 10, 50, 100, 250)
colnames(pval_table) <- sizes
rownames(pval_table) <- traits
pval_table_bb <- pval_table

print("plaque predictions... this may take over an hour")
for(k in 1:length(traits)){
  trait <- traits[k]
  for(i in 1:length(sizes)){
    size <- sizes[i]
    picrust_values <- c()
    bb_values <- c()
    for(v in 1:30){
      sam_sub_bb <- subg_bb[sample(1:nrow(subg_bb), size=size, replace=FALSE),]
      sam_sub_picrust <- subg_picrust[sample(1:nrow(subg_picrust), size=size, replace=FALSE),]
      sam_sup_bb <- suprag_bb[sample(1:nrow(suprag_bb), size=size, replace=FALSE),]
      sam_sup_picrust <- suprag_picrust[sample(1:nrow(suprag_picrust), size=size, replace=FALSE),]
      combine_picr <- rbind(sam_sub_picrust, sam_sup_picrust)
      combine_bbr <- rbind(sam_sub_bb, sam_sup_bb)
      bb_mapr <- map[rownames(map) %in% rownames(combine_bbr),]
      bb_mapr <- bb_mapr[rownames(combine_bbr),]
      bb_mapr[,map_column] <- factor(bb_mapr[,map_column])
      pic_mapr <- map[rownames(map) %in% rownames(combine_picr),]
      pic_mapr <- pic_mapr[rownames(combine_picr),]
      pic_mapr[,map_column] <- factor(pic_mapr[,map_column])
      bb_pvalr <- kruskal.test(combine_bbr[,trait] ~ bb_mapr[,map_column])$p.value
      pic_pvalr <- kruskal.test(combine_picr[,trait] ~ pic_mapr[,map_column])$p.value
      picrust_values <- c(picrust_values, pic_pvalr)
      bb_values <- c(bb_values, bb_pvalr)
    }
    bb_values[is.na(bb_values)] <- 0.99
    picrust_values[is.na(picrust_values)] <- 0.99
    pval_table_bb[trait, i] <- mean(bb_values)
    pval_table[trait, i] <- mean(picrust_values)
  }
}

write.table(pval_table_bb, 
              "data/pval_bb_plaquetest.txt", 
              quote=FALSE, 
              sep="\t", 
              col.names=NA)  
write.table(pval_table, 
              "data/pval_picrust_plaquetest.txt", 
              quote=FALSE, 
              sep="\t", 
              col.names=NA)  

print("plaque predictions written")

#Do the same for tongue and stool
pval_table_st <- data.frame(matrix(NA, nrow = 470, ncol = 5))
colnames(pval_table_st) <- sizes
rownames(pval_table_st) <- traits
pval_table_bb_st <- pval_table

print("stool-tongue predictions... this may take over an hour ")
for(k in 1:length(traits)){
  trait <- traits[k]
  for(i in 1:length(sizes)){
    size <- sizes[i]
    picrust_values <- c()
    bb_values <- c()
    for(v in 1:30){
      sam_stool_bb <- stool_bb[sample(1:nrow(stool_bb), size=size, replace=FALSE),]
      sam_stool_picrust <- stool_picrust[sample(1:nrow(stool_picrust), size=size, replace=FALSE),]
      sam_tongue_bb <- tongue_bb[sample(1:nrow(tongue_bb), size=size, replace=FALSE),]
      sam_tongue_picrust <- tongue_picrust[sample(1:nrow(tongue_picrust), size=size, replace=FALSE),]
      combine_pic <- rbind(sam_stool_picrust, sam_tongue_picrust)
      combine_bb <- rbind(sam_stool_bb, sam_tongue_bb)
      bb_map <- map[rownames(map) %in% rownames(combine_bb),]
      bb_map <- bb_map[rownames(combine_bb),]
      bb_map[,map_column] <- factor(bb_map[,map_column])
      pic_map <- map[rownames(map) %in% rownames(combine_pic),]
      pic_map <- pic_map[rownames(combine_pic),]
      pic_map[,map_column] <- factor(pic_map[,map_column])
      bb_pval <- kruskal.test(combine_bb[,trait] ~ bb_map[,map_column])$p.value
      pic_pval <- kruskal.test(combine_pic[,trait] ~ pic_map[,map_column])$p.value
      picrust_values <- c(picrust_values, pic_pval)
      bb_values <- c(bb_values, bb_pval)
    }
    bb_values[is.na(bb_values)] <- 0.99
    picrust_values[is.na(picrust_values)] <- 0.99
    pval_table_bb_st[trait, i] <- mean(bb_values)
    pval_table_st[trait, i] <- mean(picrust_values)
  }
}

write.table(pval_table_bb_st, 
              "data/pval_bb_sttest.txt", 
              quote=FALSE, 
              sep="\t", 
              col.names=NA)  
write.table(pval_table_st, 
              "data/pval_picrust_sttest.txt", 
              quote=FALSE, 
              sep="\t", 
              col.names=NA) 
print("stool tongue predictions written")


## Make Plots

library(beeswarm)
library(RColorBrewer)
library(reshape)
library(ggplot2)
library(gridExtra)

#Load the tables you wrote, (written this way to do analysis in two parts)
bb_plaque <- read.table("data/pval_bb_plaquetest.txt", sep ="\t", 
                        header=T, as.is=TRUE, check.names=FALSE)
pic_plaque <- read.table("data/pval_picrust_plaquetest.txt", sep ="\t", 
                         header=T, as.is=TRUE, check.names=FALSE)

bb_st <- read.table("data/pval_bb_sttest.txt", sep ="\t", 
                    header=T, as.is=TRUE, check.names=FALSE)
pic_st <- read.table("data/pval_picrust_sttest.txt", sep ="\t", 
                     header=T, as.is=TRUE, check.names=FALSE)

# bb_plaque[is.na(bb_plaque)] <- 0.99
# pic_plaque[is.na(pic_plaque)] <- 0.99
# bb_st[is.na(bb_st)] <- 0.99
# pic_st[is.na(pic_st)] <- 0.99

colnames(bb_plaque)[1] <- "pathway"
colnames(pic_plaque)[1] <- "pathway"
bb_plaque$approach <- "Single Cell"
pic_plaque$approach <- "Bag of Genes"
colnames(bb_st)[1] <- "pathway"
colnames(pic_st)[1] <- "pathway"
bb_st$approach <- "Single Cell"
pic_st$approach <- "Bag of Genes"

plaque <- rbind(bb_plaque, pic_plaque)
st <- rbind(bb_st, pic_st)

bb2 <- melt(plaque, id.vars=c("pathway", "approach"))
bb2$variable <- as.numeric(as.character(bb2$variable))
bb2$pathway <- as.character(bb2$pathway)
st2 <- melt(st, id.vars=c("pathway", "approach"))
st2$variable <- as.numeric(as.character(st2$variable))
st2$pathway <- as.character(st2$pathway)

cols <- brewer.pal(9, "Set1")

plaque_plot <- ggplot(bb2, aes(x=factor(variable), y=value, fill=approach)) + 
  geom_point(position=position_jitterdodge(dodge.width=0.9), size=2, 
             alpha=0.2, aes(color=approach)) +
  geom_boxplot(outlier.shape = NA, alpha=0.01, position=position_dodge(width=0.9)) +
  scale_color_manual(values=cols) +
  theme_classic() +
  theme(axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),legend.position="none") +
  xlab("Number of Samples") +
  ylab("mean p-value") +
  ggtitle("Sub- vs Supra-gingival Plaque")

st_plot <- ggplot(st2, aes(x=factor(variable), y=value, fill=approach)) + 
  geom_point(position=position_jitterdodge(dodge.width=0.9), size=2, 
             alpha=0.2, aes(color=approach)) +
  geom_boxplot(outlier.shape = NA, alpha=0.01, position=position_dodge(width=0.9)) +
  scale_color_manual(values=cols) +
  theme_classic() +
  theme(axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),legend.position="none") +
  xlab("Number of Samples") +
  ylab("mean p-value") +
  ggtitle("Stool vs Tongue")

p <- grid.arrange(st_plot, plaque_plot, ncol=2)
ggsave("Figure_3E.pdf", p, width=8, height=3)

plaque_compare <- c()
for(i in 2:(ncol(bb_plaque)-1)){
  comparison <- t.test(bb_plaque[,i], pic_plaque[,i])$p.value
  plaque_compare <- c(plaque_compare, comparison)
}
names(plaque_compare) <- c(5,10,50,100,250)

st_compare <- c()
for(i in 2:(ncol(bb_st)-1)){
  comparison <- t.test(bb_st[,i], pic_st[, i])$p.value
  st_compare <- c(st_compare, comparison)
}
names(st_compare) <- c(5,10,50,100,250)

plaque_compare2 <- c()
for(i in 2:(ncol(bb_plaque)-1)){
  comparison <- wilcox.test(bb_plaque[,i], pic_plaque[,i])$p.value
  plaque_compare2 <- c(plaque_compare2, comparison)
}
names(plaque_compare2) <- c(5,10,50,100,250)

st_compare2 <- c()
for(i in 2:(ncol(bb_st)-1)){
  comparison <- wilcox.test(bb_st[,i], pic_st[, i])$p.value
  st_compare2 <- c(st_compare2, comparison)
}
names(st_compare2) <- c(5,10,50,100,250)

sink(file="Figure_3E_Stats.txt")
cat("t-test output:\n")
cat("plaque comparisons pvals:\n")
print(plaque_compare)
cat("\n")
cat("stool-tongue comparisons pvals:\n")
print(st_compare)

cat("\nWilcox output:\n")
cat("plaque comparisons pvals:\n")
print(plaque_compare2)
cat("\n")
cat("stool-tongue comparisons pvals:\n")
print(st_compare2)

sink()
