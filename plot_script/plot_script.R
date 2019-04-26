#This script is plotting the heatmaps for output generated from motif enrichment analysis:

#output: returns heatmap plots for enriched motifs, with information of whether the gene (SF) is upregulated or downregulated

## Usage: Rscript plot_script.R  <1.INPUT_TAB_FILE> <2.OUTPUT_FOLDER> <3.MOTIF_ENRICHMENT_SCORE_MIN> <4.MOTIF_ENRICHMENT_SCORE_MAX>
## 															 <5.SF_DIFFENRETIAL_EXPRESSION_CUTOFF> <6.PLOT_FILENAME> <7.PLOT_TITLE>					

library(ggplot2)
library(grid)
library(gridExtra)
library(reshape)
library(reshape2)


##Get Arguments

args <- commandArgs(trailingOnly = TRUE) 

# ENTER ARGUMENTS HERE: (see below for description)

INPUT_FILE <- args[1]
OUTPUT_FOLDER <- args[2]
MOTIF_ENRICHMENT_SCORE_MIN <- args[3]         # >= mimimum motif enrichment score (eg. zscore)  [for plot color gradient]
MOTIF_ENRICHMENT_SCORE_MAX <- args[4]         # <= maximum zscore cutoff for plotting shade (default: 3, ie. >3 score will be rounded off to 3) 
SF_DIFFENRETIAL_EXPRESSION_CUTOFF <- args[5]  # >= minimum differntial expression score of splicing factors [for binary colour, up or down regulation, red or blue]
PLOT_FILENAME <- args[6]
PLOT_TITLE <- args[7]


#Input file format:  
#                     "COLUMNS HEADER NAME WITH DETAILS: 
#                                1.cancer	(condition type for eg. cancer type): text
#                                2.RBP_id	(Gene ID of RNA binding protein/splicing factors) : text
#                                3.regulation	 (optional field. Delta PSI value positive or negative) : binary text positive/negative
#                                4.location	 (motif presence is upstream downstream or exon) : text: up/exon/down
#                                5.robust_zscore (motif Enrichment zscore)
#                                6.exp_log2fold (RBP/SF gene Expression log2 fold change)


#### test:
#INPUT_FILE <-  "MoSEA/plot_script/test_plot.tab" #args[1] # zscore input file name 
#OUTPUT_FOLDER <- "MoSEA/plot_script" #args[2]    #"/home/babita/Desktop/MoSEA-master/test_files/plot"     # output dir to save plots

#Other parameters:
#MOTIF_ENRICHMENT_SCORE_MIN <- 1.96 	#args[3] # Minimum Zscore cutoff score (unidirectional, )
#SF_DIFFENRETIAL_EXPRESSION_CUTOFF <- 0.5 	#args[4] # log2 based differential expression score cutoff (bidirectional, absolute value)
#MOTIF_ENRICHMENT_SCORE_MAX <- 3 	#args[5] #maximum zscore for plotting shade (default: 3, ie. >3 score will be rounded off to 3) 
#PLOT_FILENAME = "test_plot.png"
#PLOT_TITLE <- "Robust-Zscore and Expression for Skipping Exons"

### ALL functions #####

#function to get the rbp_order through heatmap clustering



###################################
##Define theme
###################################

theme_heatmap <- theme_bw() + theme(plot.title=element_text(size=rel(3))) +
  theme(strip.text.x = element_text(size=16)) +
  theme(axis.ticks = element_blank()) +
  theme(axis.text.x = element_text(size = 16, angle = 45, hjust = 1, colour = "black")) +
  theme(axis.text.y = element_text(size = 16, colour = "black")) + 
  #theme(panel.grid.major.y=element_line(size=0.8),panel.grid.major.x=element_line(size=0.8)) + 
  theme(legend.position="bottom", legend.key.size=unit(1,"cm"),legend.direction='horizontal', legend.box = "horizontal", 
        legend.title = element_text(colour="black", size=12), legend.text = element_text(colour="black", size = 12)) +
  theme(plot.title = element_text(size = rel(1.5))) +
  theme(panel.spacing = unit(0.8, "cm"))  + 
  theme(plot.margin = unit(c(0, 0, 0, 0), "lines")) 


#####################################



#file without header
rzscore <-read.delim(paste0(INPUT_FILE))

#add header to the file


###2. sort dataframe in decreasing order
rzscore <- rzscore[order(rzscore$RBP_id, rzscore$robust_zscore, decreasing=TRUE),]


###3. remove duplicate rows
rzscore <- rzscore[ !duplicated(rzscore[c(1:6)]) ,]


#rzscore$regulation <- ifelse(rzscore$robust_zscore > 0, "positive", "negative")

###4.sort them back to look good
rzscore <- rzscore[order(rzscore$cancer, rzscore$RBP_id, rzscore$regulation, rzscore$location, rzscore$robust_zscore),]
######


##for cancer types & subtypes :

rzscore$cancer <- toupper(rzscore$cancer)  ##uppercase the cancer types

can_type <- unique(rzscore$cancer)

#subset data on the basis of robust_zscore > 0 and pvalue <= 0.05
rzscore_sub1 <- subset(rzscore, robust_zscore > MOTIF_ENRICHMENT_SCORE_MIN & abs(round(exp_log2fold, 2)) >= SF_DIFFENRETIAL_EXPRESSION_CUTOFF )
print(head(rzscore_sub1))
rzscore_sub1$group <- with(rzscore_sub1, ifelse(exp_log2fold >= 0,"Upregulated","Downregulated"))
rzscore_sub1$robust_zscore <- with(rzscore_sub1, ifelse(robust_zscore >= MOTIF_ENRICHMENT_SCORE_MAX, MOTIF_ENRICHMENT_SCORE_MAX, robust_zscore))     
# if z-score is greater than 3, round-off to 3. 

print(head(rzscore_sub1))
#write.table(rzscore_sub1,tmp_out,quote = FALSE,sep = "\t")


rzscore_sub1$location   <- factor(rzscore_sub1$location, levels = c("up", "exon", "down"), labels = c("Upstream", "Exon", "Downstream")) 
rzscore_sub1$regulation <- factor(rzscore_sub1$regulation, levels = c("positive", "negative"), labels = c("Positive", "Negative")) 

##write.table(rzscore_sub1,file="robust.zscore.filtered1.RBPs.table",quote=FALSE,sep="\t") #not required. #testing purposes

rzscore_sub1$robust_zscore2 <- with(rzscore_sub1, ifelse(group == "Upregulated", robust_zscore, -(robust_zscore)))  

##creating the negative z-score for heatmap clustering purpose only


#Seperate for all cancer types except brca-subtypes & exclude BLCA
##rzscore_all <- rzscore_sub1[(rzscore_sub1$cancer=="LUSC"),]# subset(rzscore_sub1, cancer!="BLCA" & (cancer == can_type_all))

#rzscore_all$cancer <- "SCLC"

rzscore_all <- rzscore_sub1

#####Run loop one for positive and one for negative subsets #for robust-zscore data only
#### 1. Plot Robust z-score heatmap on wplot

##################################
# for all cancer except subtypes
##################################

reg <- levels(rzscore_all$regulation)

i <- 1


#for(r in reg){  #use only in case of delta PSI regulation value present

######save zscores heatmap on zplot 

##### for all cancer except brca subtypes

### call function to get RBP order for each data.frame all_cancer & sub_type
#rzscore_sub <- subset(rzscore_all, regulation==r) ##use only in case of delta PSI regulation value present #to create dataframe with regulation type
rzscore_sub <- rzscore_all
rzscore_sub$RBP_id <- with(rzscore_sub, factor(RBP_id, levels=rev(levels(RBP_id))))  #re-arrange the dataframe


#reg <- ifelse(r=="Positive"," >= 0.1", " <= -0.1")



zplot <- ggplot(data=rzscore_sub, aes(x=cancer,y=RBP_id)) +  geom_tile(aes(fill=robust_zscore2), colour = "white",size=1) + 
         scale_fill_gradient2(low="midnightblue", mid="white", high="darkred", limits=c(-3,3), breaks=c(-3,0,3), 
         labels = c(paste("3 \n Downregulated"),0,paste("3 \n Upregulated") ) ) +  
         facet_grid(~location, drop = FALSE) + xlab("") + ylab("") + labs(fill = " Z-score ") + 
         theme_heatmap + ggtitle(PLOT_TITLE) #+ ggtitle(bquote(paste(Delta,psi,.(reg)))) 


#assign(paste0("zplot",i),zplot)
#print(paste0("zplot",i))


#i <- i+1

#}

#print(zplot)

#plot_path <- file.path(OUTPUT_FOLDER, PLOT_FILENAME)
ggsave(filename = PLOT_FILENAME, plot = zplot, path = OUTPUT_FOLDER, width = 7, height = 9, units = "in")
#png(plot_path,width=700, height=900)
#zplot 
graphics.off()


print("-------------")
print(paste("plot created at: ", plot_path))
#call function to create plots
#myfunction_gridPlot(zplot1, zplot2, OUTPUT_FOLDER, PLOT_FILENAME, PLOT_TITLE)

#myfunction_gridPlot(zplot, OUTPUT_FOLDER, PLOT_FILENAME, PLOT_TITLE)


# ------------------------------------ end #1 robust zscore plots -----------------------------

