
## libraries ###################################################################
require(readr)
require(ggplot2)
require(tidyr)
require(dplyr)
require(plotly)
require(WGCNA)
require(flashClust)

## functions ###################################################################

# median location normalization 
normalize.median=function (df){
  df_med = apply(as.matrix(df), 2, function(x) median(x, na.rm = T))
  df_med_loc = df_med - median(df_med)
  df.median = t(t(df)-df_med_loc)
  return(df.median)
}

## import data #################################################################

df <- read.csv("Coho DIA prot abundance Sep 2022.csv", row.names=1)
df<-t(df)
#remove fish 01b, 04, 25
df<-subset(df, select=-c(salmonovary_01b, salmonovary_04, salmonovary_25))
meta.dt <- read_csv("coho metadata Dec2022.csv")

rnames = as.matrix(rownames(df))
#df = df[,-1]
df[df == 0] = NA
df.temp = as.data.frame(df)

## missing/0 data ##############################################################

row.plot <- df.temp %>%
  mutate(id = row_number()) %>%
  gather(-id, key = "key", value = "val") %>%
  mutate(isna = is.na(val))%>%
  ggplot(aes(key, id, fill = isna)) +
  geom_raster(alpha=0.8) +
  scale_fill_manual(name = "",
                    values = c('steelblue', 'tomato3'),
                    labels = c("Present", "Missing"))+
  labs(x = "Sample",
       y = "Row Number", title = "Missing values") +
  theme(axis.text.x=element_text(angle = -65, hjust = 0, vjust = 0, size = rel(0.8)))

row.plot

## log2 transform ##############################################################

df.log = log2(df)
df.median = normalize.median(df.log)

par(mfrow = c(3,1))
boxplot(df, main = "raw")
boxplot(df.log, main = "log2 scale")
boxplot(df.median, main = "log2 scale + median normalization")

## PCA Analysis ################################################################

data.matrix = df.median

data.matrix[is.na(data.matrix)] = 0

res <- svd(data.matrix-rowMeans(data.matrix))

pca <- tryCatch({
  prcomp(t(data.matrix), retx = TRUE, center = TRUE, scale = TRUE)
}, error = function(err) {
  prcomp(t(data.matrix), retx = TRUE, center = TRUE, scale = FALSE)
})
x=1
y=2
z=3
pcVar <- round((res$d^2)/sum(res$d^2) * 100, 2)

xlab <- sprintf(paste0("PC",x,": %.2f%% var"), pcVar[x])
ylab <- sprintf(paste0("PC",y,": %.2f%% var"), pcVar[y])
zlab <- sprintf(paste0("PC",z,": %.2f%% var"), pcVar[z])

# PCA Month
pc <- pca$x
pc = cbind(as.data.frame(pc),batch=as.factor(meta.dt$Month))

fig <- plot_ly(pc, x = ~PC1, y = ~PC2, z = ~PC3, color = ~batch, 
               text= ~paste0("Month ", batch,":",rownames(pc)), hoverinfo="text")
fig <- fig %>% add_markers()
fig <- fig %>% layout(title = 'Pre-adjustment PCA (color by Month)',
                      scene = list(xaxis = list(title = xlab),
                                   yaxis = list(title = ylab),
                                   zaxis = list(title = zlab)))
fig <- fig %>% layout(legend=list(title=list(text='<b> Month </b>')))

fig
#Export -> Save as webpage

# PCA Weight
pc <- pca$x
pc = cbind(as.data.frame(pc),batch=as.factor(meta.dt$FishWeight))

fig <- plot_ly(pc, x = ~PC1, y = ~PC2, z = ~PC3, color = ~batch, 
               text= ~paste0("FishWeight", batch,":",rownames(pc)), hoverinfo="text")
fig <- fig %>% add_markers()
fig <- fig %>% layout(title = 'Pre-adjustment PCA (color by FishWeight)',
                      scene = list(xaxis = list(title = xlab),
                                   yaxis = list(title = ylab),
                                   zaxis = list(title = zlab)))
fig <- fig %>% layout(legend=list(title=list(text='<b> FishWeight </b>')))

fig

#PCA by oocyte volume
pc <- pca$x
pc = cbind(as.data.frame(pc),batch=as.factor(meta.dt$AvgVolume))

fig <- plot_ly(pc, x = ~PC1, y = ~PC2, z = ~PC3, color = ~batch, 
               text= ~paste0("AvgVolume", batch,":",rownames(pc)), hoverinfo="text")
fig <- fig %>% add_markers()
fig <- fig %>% layout(title = 'Pre-adjustment PCA (color by AvgVolume)',
                      scene = list(xaxis = list(title = xlab),
                                   yaxis = list(title = ylab),
                                   zaxis = list(title = zlab)))
fig <- fig %>% layout(legend=list(title=list(text='<b> AvgVolume </b>')))

fig


## WGCNA Prep ##################################################################

# expression matrix
datExpr = df.median
rownames(datExpr) = rnames
datExpr = as.data.frame(t(datExpr)) 

# check variance=0 and outliers
gsg = goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK #FALSE

if(!gsg$allOK)
{
    if(sum(!gsg$goodGenes)>0)
        printFlush(paste("Removing genes:", paste(names(datExpr) [!gsg$goodGenes], collapse=", ")));
    if(sum(!gsg$goodSamples)>0)
        printFlush(paste("Removing samples:", paste(rownames(datExpr) [!gsg$goodSamples], collapse=", ")));
    datExpr=datExpr[gsg$goodSamples, gsg$goodGenes]
}

# trait data
datTraits = meta.dt
rownames(datTraits) = datTraits$`File Name`

# check name alignment
rownames(datTraits)<-sub("20211005_", "", rownames(datTraits))
rownames(datTraits)<-sub(".raw", "", rownames(datTraits))
table(rownames(datTraits)==rownames(datExpr))

# fix trait data order
#datTraits = datTraits[order(rownames(datTraits)), ]
#rownames(datTraits) = datTraits$`File Name`
#table(rownames(datTraits)==rownames(datExpr)) #TRUE 61

# export data/trait
write.csv(t(datExpr), file = "datExprCohoDec2022.csv", col.names = T, row.names = T)
write.csv(datTraits, file = "datTraitsCohoDec2022.csv", col.names = T, row.names = T)

# check outlier
A = adjacency(t(datExpr),type="signed") # this calculates the whole network connectivity
k = as.numeric(apply(A,2,sum))-1 # standardized connectivity
Z.k = scale(k)
thresholdZ.k = -2.5 # often -2.5
outlierColor = ifelse(Z.k<thresholdZ.k,"red","black")
sampleTree = flashClust(as.dist(1-A), method = "average")
cnames = names(datTraits)
rnames = rownames(datTraits)
datTraits = as.data.frame(sapply(datTraits, as.character))
for (i in 1:ncol(datTraits))
  datTraits[,i] = as.numeric(as.factor(datTraits[,i]))
datTraits[is.na(datTraits == T)]=0
traitColors = data.frame(numbers2colors(datTraits,signed=FALSE))
dimnames(traitColors)[[2]] = paste(names(datTraits))
datColors = data.frame(outlier = outlierColor,traitColors)

plotDendroAndColors(sampleTree,groupLabels=names(datColors),
                    colors=datColors,main="Sample dendrogram and trait heatmap")


# remove outlying samples 
#remove.samples= Z.k<thresholdZ.k | is.na(Z.k)
#datExpr0=datExpr[!remove.samples,]
#datTraits=datTraits[!remove.samples,]
#A=adjacency(t(datExpr0),type="distance")
#k=as.numeric(apply(A,2,sum))-1
#Z.k=scale(k)
#dim(datExpr0)
#dim(datTraits)

# export data/trait
#write.csv(t(datExpr0), file = "datExprCohoNoOut.csv", col.names = T, row.names = T)
#write.csv(datTraits, file = "datTraitsCohoNoOut.csv", col.names = T, row.names = T)

# soft power threshold 
powers = c(seq(from =1, to=30, by=1)) 
sft = pickSoftThreshold(datExpr, powerVector=powers, verbose =5, 
                        networkType="signed") #signed
cex1=0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], 
     xlab= "Soft Threshold (power)", 
     ylab="Scale Free Topology Model Fit, signed R^2", type= "n", 
     main= paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], 
     labels=powers, cex=cex1, col="red")
abline(h=0.90, col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab= "Soft Threshold (power)", 
     ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")

# adjacency into topological overlap matrix
softPower = 9
adjacency = adjacency(datExpr, power = softPower, type = "signed") 
TOM = TOMsimilarity(adjacency, TOMType="signed") 
save(TOM, file="TOM.sp9CohoDec2022.RData", compress = F)


