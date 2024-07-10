library(limma)
library(tidyr)
library(ggplot2)

coho.dia<-read.csv('MSstats Input coho Sep 2022.csv')

#format the data
coho.tic2<-subset(coho.dia, select=c(File.Name, Protein.Name, Normalized.Area))
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
prot.t.tic<-prot.t.tic[2:26,]
coho.rows.tic<-rownames(prot.t.tic)
prot.t.tic<-as.data.frame(sapply(prot.t.tic, as.numeric))
rownames(prot.t.tic)<-coho.rows.tic

#check to see if normally distributed
df_tidy<-gather(prot.t.tic)

ggplot(df_tidy, aes(x=log(value))) +
  geom_density()
#data normally distributed

coho.t<-t(prot.t.tic)

#Pairwise comparison of months for differentially abundant proteins
Apr.Aug.dat<-subset(coho.t, select=c(salmonovary_10,salmonovary_11,salmonovary_13,salmonovary_16,salmonovary_22,salmonovary_28,salmonovary_17,salmonovary_20,salmonovary_24, salmonovary_26,salmonovary_27))
Aug.Dec.dat<-subset(coho.t, select=c(salmonovary_17,salmonovary_20,salmonovary_24,salmonovary_26,salmonovary_27,salmonovary_03,salmonovary_12,salmonovary_15,salmonovary_18,salmonovary_19))
Dec.Feb.dat<-subset(coho.t, select=c(salmonovary_03,salmonovary_12,salmonovary_15,salmonovary_18,salmonovary_19,salmonovary_02,salmonovary_07,salmonovary_09,salmonovary_14,salmonovary_21,salmonovary_23))

#design matrix: april vs. august
Apr.Aug.smpl<-c('salmonovary_10','salmonovary_11','salmonovary_13','salmonovary_16','salmonovary_22','salmonovary_28','salmonovary_17','salmonovary_20','salmonovary_24','salmonovary_26','salmonovary_27')
Apr.Aug.mn<-c(rep('Apr',6), rep('Aug',5))
meta.Apr.Aug<-cbind(Apr.Aug.mn, Apr.Aug.smpl)
colnames(meta.Apr.Aug)<-c('Condition', "Replicate")
meta.Apr.Aug<-data.frame(meta.Apr.Aug)

design.Apr.Aug <- model.matrix(~Condition, data = meta.Apr.Aug)
colnames(design.Apr.Aug) = gsub("Condition", "", colnames(design.Apr.Aug))

fit.Apr.Aug <- limma::lmFit(Apr.Aug.dat, design = design.Apr.Aug)

fit.Apr.Aug1 <- limma::eBayes(fit.Apr.Aug, robust = TRUE, trend = TRUE)

#find differentially abundant proteins with a p-value cut off of 0.05
diffTab.Apr.Aug = limma::topTable(fit.Apr.Aug1,
                          n=Inf, adjust = "fdr", p.value = 0.05, sort.by = "B")
write.csv(diffTab.Apr.Aug, "Limma diff abund coho proteins AprvsAug med norm.csv", quote=F)

volcanoplot(fit.Apr.Aug1, coef=2, main='April vs. August')


#design matrix: august vs. december
Aug.Dec.smpl<-c('salmonovary_17','salmonovary_20','salmonovary_24', 'salmonovary_26','salmonovary_27','salmonovary_03','salmonovary_12','salmonovary_15','salmonovary_18','salmonovary_19')
Aug.Dec.mn<-c(rep('Aug',5), rep('Dec',5))
meta.Aug.Dec<-cbind(Aug.Dec.mn, Aug.Dec.smpl)
colnames(meta.Aug.Dec)<-c('Condition', "Replicate")
meta.Aug.Dec<-data.frame(meta.Aug.Dec)

design.Aug.Dec <- model.matrix(~Condition, data = meta.Aug.Dec)
colnames(design.Aug.Dec) = gsub("Condition", "", colnames(design.Aug.Dec))

fit.Aug.Dec <- limma::lmFit(Aug.Dec.dat, design = design.Aug.Dec)

fit.Aug.Dec1 <- limma::eBayes(fit.Aug.Dec, robust = TRUE, trend = TRUE)

diffTab.Aug.Dec = limma::topTable(fit.Aug.Dec1,
                                  n=Inf, adjust = "fdr", p.value = 0.05, sort.by = "B")
write.csv(diffTab.Aug.Dec, "Limma diff abund coho proteins AugvsDec med norm.csv", quote=F)

volcanoplot(fit.Aug.Dec1, coef=2, main='August vs. December')

#design matrix: december vs. february
Dec.Feb.smpl<-c('salmonovary_03','salmonovary_12','salmonovary_15','salmonovary_18','salmonovary_19','salmonovary_02','salmonovary_07','salmonovary_09','salmonovary_14','salmonovary_21','salmonovary_23')
Dec.Feb.mn<-c(rep('Dec',5), rep('Feb',6))
meta.Dec.Feb<-cbind(Dec.Feb.mn, Dec.Feb.smpl)
colnames(meta.Dec.Feb)<-c('Condition', "Replicate")
meta.Dec.Feb<-data.frame(meta.Dec.Feb)

design.Dec.Feb <- model.matrix(~Condition, data = meta.Dec.Feb)
colnames(design.Dec.Feb) = gsub("Condition", "", colnames(design.Dec.Feb))

fit.Dec.Feb <- limma::lmFit(Dec.Feb.dat, design = design.Dec.Feb)

fit.Dec.Feb1 <- limma::eBayes(fit.Dec.Feb, robust = TRUE, trend = TRUE)

diffTab.Dec.Feb = limma::topTable(fit.Dec.Feb1,
                                  n=Inf, adjust = "fdr", p.value = 0.05, sort.by = "B")

write.csv(diffTab.Dec.Feb, "Limma diff abund coho proteins DecvsFeb med norm.csv", quote=F)

volcanoplot(fit.Dec.Feb1, coef=2, main='December vs. February')

#Compare LPN fish sampled in August.
#Three fish were outliers in the NMDS and smaller than the other 5. Find the differentially abundant proteins between these groups.
Aug.dat<-subset(coho.t, select=c(salmonovary_01b,salmonovary_04,salmonovary_17,salmonovary_20,salmonovary_24,salmonovary_25,salmonovary_26,salmonovary_27))

small.fish<-subset(coho.t, select=c(salmonovary_01b,salmonovary_04, salmonovary_25))
big.fish<-subset(coho.t, select=c(salmonovary_17,salmonovary_20,salmonovary_24,salmonovary_26,salmonovary_27))
small.sum<-rowSums(small.fish)
big.sum<-rowSums(big.fish)
sum.fish<-cbind(small.sum, big.sum)
rownames(sum.fish)<-rownames(Aug.dat)
sum.fish<-data.frame(sum.fish)

#design matrix: august
Aug.smpl<-c('salmonovary_01b','salmonovary_04','salmonovary_17','salmonovary_20','salmonovary_24','salmonovary_25', 'salmonovary_26','salmonovary_27')
weight.comp<-c('L', "L", 'H', 'H', 'H', 'L', 'H', 'H')
meta.weight<-cbind(weight.comp, Aug.smpl)
colnames(meta.weight)<-c('Condition', "Replicate")
meta.weight<-data.frame(meta.weight)

design.Aug <- model.matrix(~Condition, data = meta.weight)
colnames(design.Aug) = gsub("Condition", "", colnames(design.Aug))

fit.Aug <- limma::lmFit(Aug.dat, design = design.Aug)

fit.Aug1 <- limma::eBayes(fit.Aug, robust = TRUE, trend = TRUE)

diffTab.Aug = limma::topTable(fit.Aug1,
                              n=Inf, adjust = "fdr", p.value = 0.05, sort.by = "B")
write.csv(diffTab.Aug, "Limma diff abund coho proteins August weight.csv", quote=F)

volcanoplot(fit.Aug1, coef=2, main='small vs. big')
