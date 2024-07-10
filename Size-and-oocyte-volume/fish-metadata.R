#R code for plotting coho length, weight, condition factor, and oocyte volume by month/stage
library(ggplot2)

coho.meta<-read.csv('coho Replicates metadata2.csv')

#plot fish length by month/stage
fishlength<-ggplot(coho.meta)+
  geom_boxplot(aes(x=Month, y=FishLength, fill=Month)) +
  theme_bw() +
  scale_fill_manual(values=c("#F1EEF6", "#D7B5D8", "#DF65B0", "#CE1256"), labels=c('EPN', 'LPN', 'ECA', 'LCA'))

fishlength

#Plot fish weight by month/stage
fishweight<-ggplot(coho.meta)+
  geom_boxplot(aes(x=Month, y=FishWeight, fill=Month)) +
  theme_bw() +
  scale_fill_manual(values=c("#F1EEF6", "#D7B5D8", "#DF65B0", "#CE1256"), labels=c('EPN', 'LPN', 'ECA', 'LCA')) +
  theme(axis.text=element_text(size=15))

fishweight

#plot fish condition factor by month/stage
cond.fact<-read.csv('coho_condition_factor.csv')
cond.fact$Month<-factor(cond.fact$Month, levels=c('April', 'June', 'August'))

fishCF<-ggplot(cond.fact)+
  geom_boxplot(aes(x=Month, y=CF, fill=Month)) +
  theme_bw() +
  scale_fill_manual(values=c("#F1EEF6", "#D7B5D8", "#DF65B0"), labels=c('EPN', 'LPN', 'ECA'))

fishCF

#plot oocyte volume by month/stage
oocyte.vol<-read.csv('coho_egg_volumes.csv')
oocyte.vol$Month<-factor(oocyte.vol$Month, levels=c('April', 'June', 'August', 'October', 'December', 'February'))

fishvol<-ggplot(oocyte.vol)+
  geom_boxplot(aes(x=Month, y=Volume, fill=Month)) +
  theme_bw() +
  scale_fill_manual(values=c("#F1EEF6", "#D7B5D8", "#DF65B0", "#CE1256"), labels=c('EPN', 'LPN', 'ECA', 'LCA'))+
  theme(axis.text=element_text(size=13))

fishvol
