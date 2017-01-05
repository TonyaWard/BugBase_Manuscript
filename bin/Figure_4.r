#### Create Files Needed for Figure 4 ####
library(RColorBrewer)

#Make a list of the stool OTU IDs to filter the green genes tree by
contributions <- read.table("Stool_only/otu_contributions/contributing_otus.txt", 
                  sep="\t", 
                  comment="",
                  row=1,
                  header=T, 
                  as.is=T, 
                  check.names=F)
#Convert from logic to integer
contributions <- contributions * 1


taxa <- read.table("data/97_otu_taxonomy.txt.gz", 
                   sep="\t", 
                   comment="",
                   row=1,
                   as.is=T, 
                   check.names=F)
taxa <- taxa[rownames(contributions),,drop=F]

#split taxonomy at different levels
otu_names <- as.character(taxa$V2)
names(otu_names) <- rownames(taxa)
names <- read.table(text = otu_names, sep = ";", colClasses = "character")
rownames(names) <- names(otu_names)

#Remove the level annotations
names$V1 <- gsub('k__', '', names$V1)
names$V2 <- gsub(' p__', '', names$V2)
names$V3 <- gsub(' c__', '', names$V3)
names$V4 <- gsub(' o__', '', names$V4)
names$V5 <- gsub(' f__', '', names$V5)
names$V6 <- gsub(' g__', '', names$V6)
names$V7 <- gsub(' s__', '', names$V7)
colnames(names) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

#Keep just one OTU per genus
names <- names[!duplicated(names$Genus),]

#Subset contributions to these OTU IDs
contributions <- contributions[rownames(names),]
#Write to table
write.table(contributions, file="data/Stool_Input.txt", quote=F, sep="\t", col.names=NA)

#Write these IDs to a file
write(rownames(contributions), file = "data/Stool_IDs.txt", ncolumns=1)


#Write in a color for each unique phylum
phyla <- unique(names$Phylum)
#Get 38 unique colors
cols <- brewer.pal(9,'Set1')
cols2 <- colorRampPalette(cols)(length(phyla))

for(i in 1:length(phyla)){
  phylum <- phyla[i]
  for(r in 1:nrow(names)){
    if(names[r,"Phylum"] == phylum){
      names[r,"Species"] <- cols2[i]
    }
  }
}

#Write to file as stool_taxa_map.txt
write.table(names, file="data/stool_taxa_map.txt", quote=F, sep='\t', col.names=NA)

print("Writing annotation file...")

#Write the graphlan annotation file
sink("data/annotation.txt")
cat("title	HMP Stool Bacteria
title_font_size 20
total_plotted_degrees	340
annotation_background_alpha	0.1
start_rotation	270
branch_bracket_width	0.0
class_legend_font_size	14
annotation_legend_font_size	24	
annotation_background_separation	-0.03
annotation_background_offset	0.0
annotation_background_width	0.07
branch_thickness	1.5
ring_label_font_size	1	10
ring_label_font_size	2	10
ring_label_font_size	3	10
ring_label_font_size	4	10
ring_label_font_size	5	10
ring_label_font_size	6	10
clade_marker_size	10
clade_marker_edge_width	0.8
ring_internal_separator_thickness	1	0.5
ring_internal_separator_thickness	2	0.5
ring_internal_separator_thickness	3	0.5
ring_internal_separator_thickness	4	0.5
ring_internal_separator_thickness	5	0.5
ring_internal_separator_thickness	6	0.5
ring_label_color	1	black
ring_label_color	2	black
ring_label_color	3	black
ring_label_color	4	black
ring_label_color	5	black
ring_label_color	6	black")
sink()

#Add annotations for each OTU
for(r in 1:nrow(contributions)){
  row <- rownames(contributions)[r]
  sink("data/annotation.txt", append=T)
  #set gram staining (black = positive)
  cat(paste("\n",row,"\tring_color\t1\tblack",sep=""))
  if(contributions[row,"Gram_Negative"] > 0){
    cat(paste("\n",row,"\tring_alpha\t1\t0.2",sep=""))
  }
  #Set oxygen tolerance colors (blue = aerobic)
  cat(paste("\n",row,"\tring_color\t2\tblue",sep=""))
  if(contributions[row,"Facultatively_Anaerobic"] > 0){
    cat(paste("\n",row,"\tring_alpha\t2\t0.25",sep=""))
  }
  if(contributions[row,"Anaerobic"] > 0){
    cat(paste("\n",row,"\tring_alpha\t2\t0.01",sep=""))
  }
  #set stress tolerance
  cat(paste("\n",row,"\tring_color\t3\torange",sep=""))
  if(contributions[row,"Stress_Tolerant"] == 0){
    cat(paste("\n",row,"\tring_alpha\t3\t0.2",sep=""))
  }
  #set mobile elements
  cat(paste("\n",row,"\tring_color\t4\tgreen",sep=""))
  if(contributions[row,"Contains_Mobile_Elements"] == 0){
    cat(paste("\n",row,"\tring_alpha\t4\t0.2",sep=""))
  }
  #set biofilm formation
  cat(paste("\n",row,"\tring_color\t5\tyellow",sep=""))
  if(contributions[row,"Forms_Biofilms"] == 0){
    cat(paste("\n",row,"\tring_alpha\t5\t0.2",sep=""))
  }
  #set pathogenicity
  cat(paste("\n",row,"\tring_color\t6\tred",sep=""))
  if(contributions[row,"Pathogenic"] == 0){
    cat(paste("\n",row,"\tring_alpha\t6\t0.2",sep=""))
  }
  sink()
}

for(i in 1:nrow(names)){
  row <- rownames(names)[i]
  sink("data/annotation.txt", append=T)
  cat(paste("\n",row,"\tclade_marker_color\t", names[row,ncol(names)],sep=""))
  cat(paste("\n",names[row,"Phylum"], "\tannotation\t", names[row,"Phylum"],sep=""))
  cat(paste("\n",names[row,"Phylum"], "\tclade_marker_color\t", names[row,ncol(names)],sep=""))
  sink()
}

print("Annotation file complete.")
