library(vegan)

#matrix_file<-"~/Repos/CBAS/CBAS_16S/Temperature_RawFiles/cbas_tempVSctrl_NoOTUCountLess50.otutab.csv"
#matrix_file<-"~/Repos/CBAS/CBAS_16S/Temperature_RawFiles/cbas_tempVSctrl.otutab.bactloadCorrected_valuesOnly_wFreqs_more50Counts.csv"
matrix_file<-"~/Repos/CBAS/CBAS_16S/Temperature_RawFiles/cbas_tempVSctrl.otutab.bactloadCorrected_valuesOnly_wFreqs_more1PctCounts.csv"

bact_matrix<-read.csv(matrix_file, head=T, sep="\t", row.names=1)
colnames(bact_matrix)

#groups = GreenControl = 1
#         GreenTemp = 2
#         PurpleControl = 3
#         PurpleTemp = 4
#num_groups<-c(1,1,1,2,2,1,1,1,2,3,3,3,4,4,4,3,3,3,4,4)
#label_groups<-c("GC","GC","GC","GT","GT","GC","GC","GC","GT","PC","PC","PC","PT","PT","PT","PC","PC","PC","PT","PT")

#wFreqs_more50Counts % more1PctCounts
num_groups<-c(1,1,1,2,2,1,1,1,2,3,3,4,3,3,3,4)
label_groups<-c("GC","GC","GC","GT","GT","GC","GC","GC","GT","PC","PC","PT","PC","PC","PC","PT")

#head(bact_matrix)

nmds<-metaMDS(t(bact_matrix), try=50)
plot(nmds, type = "t", display="sites")

points(nmds, col=num_groups, cex=1.5, pch=16)

ordispider(nmds, label_groups, label=T)

anosim(t(bact_matrix), num_groups)


#temp_treatments<-data.frame(Condition=c("C","C","C","T","T","C","C","C","T","C","C","C","T","T","T","C","C","C","T","T")
#                            , Sponge_type=c("G","G","G","G","G","G","G","G","G","P","P","P","P","P","P","P","P","P","P","P"))

#wFreqs_more50Counts % more1PctCounts
temp_treatments<-data.frame(Condition=c("C","C","C","T","T","C","C","C","T","C","C","T","C","C","C","T")
                            , Sponge_type=c("G","G","G","G","G","G","G","G","G","P","P","P","P","P","P","P"))

adonis(t(bact_matrix)~Condition, data=temp_treatments, strata=temp_treatments$Sponge_type)

bact_cca<-cca(t(bact_matrix)~Condition+Sponge_type+Condition:Sponge_type, data=temp_treatments)

plot(bact_cca, type="t", display="sites")
points(bact_cca, col=num_groups, cex=1.5, pch=16)
ordispider(bact_cca)
anova(bact_cca)
anova(bact_cca, by="term")


