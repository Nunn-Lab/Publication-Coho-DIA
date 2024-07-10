library(tidyr)
library(vegan)
library(WGCNA)
source('biostats.R')
library(RColorBrewer)

#NMDS of coho DIA proteins

#This is a large file and will take a moment to read into R
coho.dia<-read.csv('MSstats Input coho Sep 2022.csv')

area.med<-median(coho.dia$Area)

######Workflow for the entire dataset. See below for the workflow after removal of outlier fish.
#select color palette
brewer.pal(n=4, name='PuRd')

#create grouping by dominant ovarian follicle stage for each time point
group.coho<-c('LPN', 'LCA', 'ECA', 'LPN', 'LCA','LCA','EPN','EPN','ECA','EPN','LCA','ECA','EPN','LPN','ECA','ECA','LPN','LCA','EPN','LCA','LPN','LPN','LPN','LPN','EPN')

#create data frame for plotting
coho.tic2<-subset(coho.dia, select=c(File.Name, Protein.Name, Normalized.Area))
coho.tic2$Normalized.Area<-as.numeric(coho.tic2$Normalized.Area)

#sum peptide peak areas to create aggregate peak area for each protein
ProtSums.tic<-aggregate(. ~ File.Name+Protein.Name, data=coho.tic2, FUN=sum)
#reformat so that column File.Name is column headers and Protein.Name is row names
ProtSums.wide.tic<-pivot_wider(ProtSums.tic, names_from=File.Name, values_from=Normalized.Area)

ProtSums.df.tic<-data.frame(ProtSums.wide.tic)

rownames(ProtSums.df.tic)<-ProtSums.df.tic$Protein.Name
colnames(ProtSums.df.tic)<-sub('X20211005_', '', colnames(ProtSums.df.tic))
colnames(ProtSums.df.tic)<-sub('.raw', '', colnames(ProtSums.df.tic))

prot.t.tic<-as.data.frame(t(ProtSums.df.tic))
prot.t.tic<-prot.t.tic[2:26,]
coho.rows.tic<-rownames(prot.t.tic)
prot.t.tic<-as.data.frame(sapply(prot.t.tic, as.numeric))
rownames(prot.t.tic)<-coho.rows.tic

#check for outliers
gsg = goodSamplesGenes(prot.t.tic, verbose=3)
gsg$allOK

#NMDS
#log transform
coho.tra.tic<-data.trans((prot.t.tic+1), method='log', plot=F)

nmds.tic<-metaMDS(coho.tra.tic, distance='bray', trymax=10, autotransform=F)

#creat plot
fig.tic<-ordiplot(nmds.tic, type='none', display='sites', xlab='NMDS1', ylab='NMDS2', choices=c(1,2))
points(fig.tic, 'sites', pch=19, col=c("#D7B5D8", "#CE1256", "#DF65B0", "#D7B5D8", rep("#CE1256",2),rep("#F1EEF6",2),"#DF65B0","#F1EEF6","#CE1256","#DF65B0","#F1EEF6","#D7B5D8", rep("#DF65B0",2),"#D7B5D8","#CE1256","#F1EEF6","#CE1256", rep("#D7B5D8",4),"#F1EEF6"))
ordihull(fig.tic,groups=group.coho,draw='lines',col='grey75',label=T)

######Three fish sampled in August (LPN) are outliers and are smaller than the rest. Remove these from the analysis.
coho.sub<-subset(coho.dia, File.Name!='20211005_salmonovary_01b.raw')
coho.sub<-subset(coho.sub, File.Name!='20211005_salmonovary_04.raw')
coho.sub<-subset(coho.sub, File.Name!='20211005_salmonovary_25.raw')

coho.tic2<-subset(coho.sub, select=c(File.Name, Protein.Name, Normalized.Area))
coho.tic2$Normalized.Area<-as.numeric(coho.tic2$Normalized.Area)

coho.tic2$med.norm<-log((coho.tic2$Normalized.Area*3972043)+1)
coho.tic2<-subset(coho.tic2, select=c(File.Name, Protein.Name,med.norm))

ProtSums.tic<-aggregate(. ~ File.Name+Protein.Name, data=coho.tic2, FUN=sum)
#reformat so that column File.Name is column headers and Protein.Name is row names
ProtSums.wide.tic<-pivot_wider(ProtSums.tic, names_from=File.Name, values_from=med.norm)

ProtSums.df.tic<-data.frame(ProtSums.wide.tic)

rownames(ProtSums.df.tic)<-ProtSums.df.tic$Protein.Name
colnames(ProtSums.df.tic)<-sub('X20211005_', '', colnames(ProtSums.df.tic))
colnames(ProtSums.df.tic)<-sub('.raw', '', colnames(ProtSums.df.tic))

prot.t.tic<-as.data.frame(t(ProtSums.df.tic))
prot.t.tic<-prot.t.tic[2:23,]
coho.rows.tic<-rownames(prot.t.tic)
prot.t.tic<-as.data.frame(sapply(prot.t.tic, as.numeric))
rownames(prot.t.tic)<-coho.rows.tic

gsg = goodSamplesGenes(prot.t.tic, verbose=3)
gsg$allOK

#NMDS
group.coho<-c('LCA', 'ECA', 'LCA','LCA','EPN','EPN','ECA','EPN','LCA','ECA','EPN','LPN','ECA','ECA','LPN','LCA','EPN','LCA','LPN','LPN','LPN','EPN')

coho.tra.tic<-data.trans((prot.t.tic+1), method='log', plot=F)
nmds.tic<-metaMDS(coho.tra.tic, distance='bray', trymax=10, autotransform=F)

fig.tic<-ordiplot(nmds.tic, type='none', display='sites', xlab='NMDS1', ylab='NMDS2', choices=c(1,2))
points(fig.tic, 'sites', pch=19, col=c( "#CE1256", "#DF65B0",  rep("#CE1256",2),rep("#F1EEF6",2),"#DF65B0","#F1EEF6","#CE1256","#DF65B0","#F1EEF6","#D7B5D8", rep("#DF65B0",2),"#D7B5D8","#CE1256","#F1EEF6","#CE1256", rep("#D7B5D8",3),"#F1EEF6"),cex=1.5)
ordihull(fig.tic,groups=group.coho,draw='lines',col='grey30',label=T)

#####ANOSIM to determine significant of groupings
coho.row<-data.stand(prot.t.tic, method='total', margin='row', plot=F)
coho.d<-vegdist(coho.row, 'bray')
coho.mon<-anosim(coho.d, grouping=group.coho)
summary(coho.mon)
#ANOSIM statistic R: 0.8276 
#Significance: 0.001 
