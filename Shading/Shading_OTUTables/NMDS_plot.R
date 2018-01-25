library(vegan)

bact_matrix<-read.csv("~/Repos/CBAS/CBAS_16S/Shading_OTUTables/cbas_bleachVSctrl_NoLBSponges_NoOTUCountLess50.otutab.standardizedsamplenames.csv", head=T, sep="\t", row.names=1)
num_groups<-c(1,2,2,2,2,1,1,1,1,1,1,3,4,5,3,2,4,5,1)
label_groups<-c("C","B","B","B","B","C","C","C","C","C","C","W9","W3","W6","W9","B","W3","W6","C")

#head(bact_matrix)

nmds<-metaMDS(t(bact_matrix), try=50)
plot(nmds, type = "t", display="sites")

points(nmds, col=num_groups, cex=1.5, pch=16)

ordispider(nmds, label_groups, label=T)

anosim(t(bact_matrix), num_groups)
