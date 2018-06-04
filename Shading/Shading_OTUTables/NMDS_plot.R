library(vegan)
library(ggplot2)
library(alluvial)
library(reshape2)

OTU50_bact_matrix<-read.csv("~/Repos/CBAS/CBAS_16S/Shading_RawFiles/cbas_bleachVSctrl_NoLBSponges_NoOTUCountLess50.otutab.newsamplenames_bactLoadStandardized.csv", head=T, sep=",", row.names=1)
num_groups<-c(1,2,2,2,2,1,1,1,1,1,3,4,5,3,2,4,5)
label_groups<-c("C","B","B","B","B","C","C","C","C","C","W9","W3","W6","W9","B","W3","W6")

OTU1Pct_bact_matrix<-read.csv("~/Repos/CBAS/CBAS_16S/Shading_RawFiles/cbas_bleachVSctrl_NoLBSponges_NoOTULessOnePercentCNs.otutab.newsamplenames_bactLoadStandardized.csv", head=T, sep=",", row.names=1)

bact_matrix<-OTU1Pct_bact_matrix
#bact_matrix<-OTU50_bact_matrix

#head(bact_matrix)

nmds<-metaMDS(t(bact_matrix), try=50)
plot(nmds, type = "t", display="sites")

points(nmds, col=num_groups, cex=1.5, pch=16)

ordispider(nmds, label_groups, label=T)

anosim(t(bact_matrix), num_groups)

bact_load<-read.csv("~/Repos/CBAS/CBAS_16S/Shading_RawFiles/cbas_bleachVSctrl_Bact_Load.csv")

ggplot(bact_load, aes(Sample, log(V4_Copy_number))) + geom_boxplot() + geom_jitter() + theme_bw() 


##alluvial plot

OTU1Pct_bact_matrix_melted<-melt(t(OTU1Pct_bact_matrix[,c(1:10,15)]), varnames = c("Sample", "OTU"))
text_treat<-strsplit(as.character(OTU1Pct_bact_matrix_melted$Sample), "_")
treatment<-lapply(text_treat, function(x) x[2])
OTU1Pct_bact_matrix_melted<-cbind(OTU1Pct_bact_matrix_melted, as.character(treatment))
colnames(OTU1Pct_bact_matrix_melted)<-c("Sample","OTU","value","Treatment")

OTU1Pct_bact_matrix_aggregated<-aggregate(value~Treatment+OTU, data=OTU1Pct_bact_matrix_melted, sum)
OTU1Pct_bact_matrix_aggregated_ctrl<-subset(OTU1Pct_bact_matrix_aggregated, Treatment=="Control")
OTU1Pct_bact_matrix_aggregated_bleached<-subset(OTU1Pct_bact_matrix_aggregated, Treatment=="Bleached")

OTU1Pct_bact_matrix_aggregated_ctrl<-cbind(OTU1Pct_bact_matrix_aggregated_ctrl, OTU1Pct_bact_matrix_aggregated_ctrl$value/sum(OTU1Pct_bact_matrix_aggregated_ctrl$value))
colnames(OTU1Pct_bact_matrix_aggregated_ctrl)<-c("Treatment","OTU","value","Freq")
OTU1Pct_bact_matrix_aggregated_bleached<-cbind(OTU1Pct_bact_matrix_aggregated_bleached, OTU1Pct_bact_matrix_aggregated_bleached$value/sum(OTU1Pct_bact_matrix_aggregated_bleached$value))
colnames(OTU1Pct_bact_matrix_aggregated_bleached)<-c("Treatment","OTU","value","Freq")

OTU1Pct_bact_matrix_aggregated<-rbind(OTU1Pct_bact_matrix_aggregated_ctrl,OTU1Pct_bact_matrix_aggregated_bleached)

alluvial(OTU1Pct_bact_matrix_aggregated[,c(2,1)], freq=OTU1Pct_bact_matrix_aggregated$Freq, col=ifelse(OTU1Pct_bact_matrix_aggregated$Treatment == "Control",  "lightslateblue", "slategray1"), layer=OTU1Pct_bact_matrix_aggregated$Treatment=="Bleached", alpha=0.6)



