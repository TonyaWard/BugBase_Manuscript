"otu.plots.fig2" <- function(input_traits, input_otutable, input_otuscontributing, input_map, input_taxonomy){
  
	#define inputs
	traits <- input_traits
	otu_table <- input_otutable
	otus_contributing <- input_otuscontributing
	map <- input_map
	taxonomy <- input_taxonomy
	map_column <- "HMPBODYSUBSITE"

	#Convert to RA
	otu_table <- sweep(otu_table, 2, colSums(otu_table), FUN='/')
	
	#set taxonomy level
	taxa_level <- 2
	taxa_level <- as.numeric(taxa_level)

	#read in gg_taxonomy table
	gg_taxonomy <- read.table(taxonomy, sep="\t", row.names=1, check=F, quote='')
	  
	#keep only otus in the gg taxonomy table that are in the otu table
	#these tables will be in the same otu order
	#OTUs are now rows
	otus_keep <- intersect(rownames(otu_table),rownames(gg_taxonomy))
	if(length(otus_keep) < nrow(otu_table)){
		print("Some OTUs did not have a taxonomy associated. Did you pick against Greengenes 97% or IMG?")
	}
	gg_taxonomy <- gg_taxonomy[otus_keep,, drop=F]
	otu_table <- otu_table[otus_keep,,drop=F]
	otus_contributing <- otus_contributing[otus_keep,,drop=F]

	#subset the gg taxonomy to the level specified (default=2)
	names_split <- array(dim=c(length(gg_taxonomy[,1]), 7))
	otu_names <- as.character(gg_taxonomy[,1])
	for(i in 1:length(otu_names)){
		names_split[i,] <- strsplit(otu_names[i], ";", fixed=T)[[1]]
	}
	otu_names <- names_split[,taxa_level]
	for(i in 1:length(otu_names)){
		otu_names[i] <- strsplit(otu_names[i], "__", fixed=T)[[1]][2]
	}
	names_split[,taxa_level] <- otu_names
	for(i in 1:nrow(names_split)){
		if(is.na(names_split[i,taxa_level])){
			if(taxa_level > 1){
				names_split[i, taxa_level] <- names_split[i, taxa_level -1]
			} else {
				names_split[i, taxa_level] <- "unknown"
			}
		}
	}

	#store the full otu_table
	otu_table1 <- otu_table
	
	#add taxonomy as the rownames in the otu table
  otu_table <- as.matrix(otu_table)
	rownames(otu_table) <- names_split[,taxa_level]

	#aggregate to the same taxa
	#you must t() again, to have samples as columns
	otu_table <- t(sapply(by(otu_table,rownames(otu_table),colSums),identity))

	#ensure same order of samples in map and otu table
	map <- map[colnames(otu_table),,drop=F]

	#create as many colors as there are taxa, and name them with the taxa names
	#keep the colors similar to those in the other plots (for manuscript aes)
	#cols_otus <- colorRampPalette(brewer.pal(9,'Set1'))
	#cols2 <- cols_otus(length(unique(rownames(otu_table))))
	cols2 <- c("#E41A1C", "#BE2F3D", "#E41A1C", "#735A81", "#4E70A2", "#377EB8", "#3D8C97", "#42977F", "#47A167", "#4BAC4F",
	           "#599F58", "#698A6B", "#79757E", "#896092", "#9A4F9E", "#4DAF4A", "#984EA3", "#DD6F34", "#F37911", "#FF8C05",
	           "#FFA810", "#FFC41B", "#FFDF26", "#FFFB31", "#EEDF30", "#FF7F00", "#C7952C", "#B47129", "#AA5830", "#BB6150",
	           "#CD6A71", "#DE7492", "#F07DB2", "#EA84B9", "#D589B1", "#C18EA9", "#AD93A1", "#999999")
	names(cols2) <- unique(rownames(otu_table))
	cols2 <- c(cols2,"#C0C0C0")
	names(cols2)[length(cols2)] <- "Other"

	taxa_list <- c()
	#create taxa summaries
	taxa_plots <- c()
	trait_names <- gsub("_", " ", traits)
	for(x in 1:length(traits)){
	  trait <- traits[x]
	  
	  #multiple the full otu_table by the boolean table for otus contributing
	  #   to the trait
	  positive_otu_table <- t(sweep(t(otu_table1), 2, 
	                                otus_contributing[,trait],"*"))
	  
	  #add taxonomy as the rownames in the otu table
	  rownames(positive_otu_table) <- names_split[,taxa_level]
	  
	  #aggregate to the same taxa
	  #you must t() again, to have samples as columns
	  positive_otu_table <- t(sapply(by(positive_otu_table,rownames(positive_otu_table),colSums),identity))
	  
	  #ensure same order of samples in map and otu table
	  map <- map[colnames(positive_otu_table),,drop=F]
	  
	  #melt the otu_table by sample id
	  melted_otu_table <- melt(positive_otu_table)
	  colnames(melted_otu_table) <- c("Taxa", "SampleID", "Count")
	  
	  #merge the otu table and mapping file
	  map$SampleID <- rownames(map)
	  melted_otu_table <- merge(melted_otu_table, map, by="SampleID")
	  
	  #collapse by groups in the map column
	  group_collapsed_otus <- ddply(melted_otu_table, .(Taxa,melted_otu_table[,map_column]), summarize, Count = mean(Count))
	  
	  colnames(group_collapsed_otus)[2] <- map_column
	  
	  #set value for cutoff (1/10 of the highest proportion)
	  max_abund <- max(group_collapsed_otus$Count)
	  cutoff_val <- max_abund / 10
	  
	  #call taxa that are less than cutoff of the population "Other"
	  group_collapsed_otus$Taxa <- as.character(group_collapsed_otus$Taxa)
	  group_collapsed_otus$Count <- as.numeric(group_collapsed_otus$Count)
	  
	  group_collapsed_otus[which(group_collapsed_otus$Count < cutoff_val),"Taxa"] <- "Other"
	  #re-collapse to group the 'Others'
	  group_collapsed_otus <- ddply(group_collapsed_otus, 
	                                .(Taxa,group_collapsed_otus[,2]), 
	                                summarize, Count = sum(Count))
	  colnames(group_collapsed_otus)[2] <- map_column
	  
	  taxa_list <- c(taxa_list, unique(group_collapsed_otus$Taxa))
	  
    group_collapsed_otus[,map_column] <- as.character(group_collapsed_otus[,map_column])
	  group_collapsed_otus <- group_collapsed_otus[order(group_collapsed_otus[,map_column]),]
	  
	  colnames(group_collapsed_otus)[2] <- "this_column"
	  #make the plot
	  taxa_plot <- NULL
	  taxa_plot <- ggplot(group_collapsed_otus, aes_string(x = "this_column",y = "Count", fill="Taxa")) + 
	    geom_bar(stat="identity", show_guide=FALSE) + 
	    labs(y = "Relative Abundance", x = "", title=trait_names[x]) +
	    theme_classic() +
	    theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'), axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + 
	    scale_fill_manual(values=cols2)
	  
		taxa_plots[[trait]] <- taxa_plot 
	}

	#subset to the taxa that met the cutoff
	taxa_list <- unique(taxa_list)
	cols_keep <- cols2[which(names(cols2) %in% taxa_list)]

	#create a legend table with names, colors and coordinates
	legend_info <- matrix(0, length(cols_keep), 4)
	legend_info[,1] <- names(cols_keep)
	legend_info[,4] <- 1
	counter <- c(1:length(cols_keep)+1)
	for(i in 1:length(cols_keep)){
	  legend_info[i,2] <- cols_keep[[i]]
	  legend_info[i,3] <- counter[i]
	}
	colnames(legend_info) <- c("Taxa", "Color", "Y", "X")
	legend_info <- as.data.frame(legend_info)
	legend_info$X <- as.numeric(legend_info$X)
	legend_info$Y <- as.numeric(legend_info$Y)
	legend_info$Color <- as.character(legend_info$Color)
	#Add a space before name so plot legend looks nice
	legend_info$Taxa <- sub("^", " ", legend_info$Taxa ) 

	#name pdf
	taxa_legend <- c("taxa_legend_fig2.pdf")
	
	#make the pdf
	pdf(taxa_legend, height=6,width=6)
	par(mar=c(0.5,0.5,0.5,0.5), oma=c(0.1,0.1,0.1,0.1), mgp=c(1.5,0.5,0))
	
	#plot the points
	plot(legend_info$X, legend_info$Y, 
		 pch=15, 
		 cex=3, 
		 col= legend_info$Color, 
		 axes=FALSE, 
		 xlab='', 
		 ylab='')
	#add names
	text(legend_info$X, legend_info$Y, legend_info$Taxa, pos=4)
	graphics.off()
	return(taxa_plots)
}
