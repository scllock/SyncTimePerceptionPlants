

#libraries
library(ggplot2) 
library(maps) 
library(pheatmap)
library(seqinr)
library(RColorBrewer)

#################################
#  In this section, we will explore two KH genes in a range of ecotypes
#################################

kh17="AT3G32940"
kh29="AT5G56140"

############################
#Q1: how do the sequences vary over ecotypes?

#first, we downloaded the pseudogenes of these genes from the 1001 genome project
#then, we used clustalo to go a multiple sequence alignment

#Step 1: load in multiple sequence alignment
kh17Aln=read.alignment('data/ecotype/clustal0-kh17.aln-clustal_num', format="clustal")
kh29Aln=read.alignment('data/ecotype/clustal0-kh29.aln-clustal_num', format="clustal")


png('plots/Ecotypes/FigS_mutationalProfilekh1729.png', width=4, height=6, unit='in', res=600)
par(mfrow=c(2,1))


kh17SeqAln=sapply(kh17Aln$seq, function(i){
  seq=i[1]
  seq=gsub("\t", "", seq)
  seq=gsub('[[:digit:]]+', '', seq)
  seq
})

kh29SeqAln=sapply(kh29Aln$seq, function(i){
  seq=i[1]
  seq=gsub("\t", "", seq)
  seq=gsub('[[:digit:]]+', '', seq)
  seq
})

#Step 2: How much does each sequence vary from the consensus sequence?

#a) calculate consensus sequence
matChar=sapply(strsplit(kh17SeqAln, ""), function(i){unlist(i)})
kh17MatChar=matChar

freq=apply(matChar, 1, function(i){
  table(factor(i, levels=c("-", "n", "a", "c", "g", "t")))
})


consensus=apply(freq, 2, function(i){
  names(i)[which.max(i)]
})

#print sequence without gaps
paste(gsub("-", "",consensus), collapse='')

#print sequence with gaps
paste(consensus, collapse='')

#this was searched in blastx to identify the locations of conserved domains:
#KH: 1768-1857
#Atrophin-1:230-967

#calculate which regions have the most mutations
consensusLevel=apply(freq, 2, function(i){
  c(max(i)/sum(i), max(i[3:6])/sum(i[3:6]))
})

#average over space for greater clarity:
by=10
movingAvg=sapply(c(1:(length(consensusLevel[1,])-by)), function(i){
  mean(consensusLevel[1,i:(i+by)])
})
plot(movingAvg, type='l', xlab='position (bp)', ylab='identity', main="KH17 mutational profile")
#KH: 1768-1857
#Atrophin-1:230-967
abline(v=c(1768, 1857, 230, 967)+by/2)
hotspots=which(consensusLevel[2,]<0.95)
hotspots=hotspots[which(hotspots<length(movingAvg))]
points(c(1:length(movingAvg))[hotspots], consensusLevel[2,hotspots], pch=20, col='grey')


###########
##repeat for kn29
matChar=sapply(strsplit(kh29SeqAln, ""), function(i){unlist(i)})
k29MatchChar=matChar

freq=apply(matChar, 1, function(i){
  table(factor(i, levels=c("-", "n", "a", "c", "g", "t")))
})


consensus=apply(freq, 2, function(i){
  names(i)[which.max(i)]
})

#print sequence without gaps
paste(gsub("-", "",consensus), collapse='')

#print sequence with gaps
paste(consensus, collapse='')

#this was searched in blastx to identify the locations of conserved domains:
#KH 2191-2331	0e+00
#KH 1871-1915	0e+00
#KH 1992-2108	0e+00

#calculate which regions have the most mutations
consensusLevel=apply(freq, 2, function(i){
  c(max(i)/sum(i), max(i[3:6])/sum(i[3:6]))
})

#average over space for greater clarity:
by=10
movingAvg=sapply(c(1:(length(consensusLevel[1,])-by)), function(i){
  mean(consensusLevel[1,i:(i+by)])
})
plot(movingAvg, type='l', xlab='position (bp)', ylab='identity', main="KH29 mutational profile")

#KH 2191-2331	0e+00
#KH 1871-1915	0e+00
#KH 1992-2108	0e+00
abline(v=c(2191, 2331, 1871, 1915, 1992, 2108)+by/2)
hotspots=which(consensusLevel[2,]<0.95)
hotspots=hotspots[which(hotspots<length(movingAvg))]
points(c(1:length(movingAvg))[hotspots], consensusLevel[2,hotspots], pch=20, col='grey')

dev.off()



#############################################
#Q2: How are the ecotypes grouped and how does that correlate with the world map?

#Step 1: cluster mutational profiles of ecotypes

matChar=kh17MatChar
numberMatches=apply(matChar,2, function(i){
  apply(matChar,2, function(j){
    length(which(i==j & i!='n' & j!='n'))
  })
})

#include only sequences without too many ns: THRESHOLD=90%*********
notNCount=apply(matChar, 2, function(i){
  length(which(i!='n'))
})
hist(notNCount)

ids=which(notNCount>(0.9*max(notNCount)))

a=pheatmap(numberMatches[ids, ids]/dim(matChar)[1], cutree_rows=5, cutree_cols = 5, file='plots/Ecotypes/kh17_ecotypeClusters.png')

kh17Clusts=cutree(a$tree_row, 5)

kh17SubNames=kh17Aln$nam[ids]
names(kh17Clusts)=kh17SubNames

##repeat with kh29

matChar=k29MatchChar
numberMatches=apply(matChar,2, function(i){
  apply(matChar,2, function(j){
    length(which(i==j & i!='n' & j!='n'))
  })
})

#include only sequences without too many ns: THRESHOLD=90%*********
notNCount=apply(matChar, 2, function(i){
  length(which(i!='n'))
})
hist(notNCount)

ids=which(notNCount>(0.9*max(notNCount)))
a=pheatmap(numberMatches[ids, ids]/dim(matChar)[1])
a=pheatmap(numberMatches[ids, ids]/dim(matChar)[1], cutree_rows=5, cutree_cols = 5, file='plots/Ecotypes/KH29_ecotypeClusters.png')

kh29Clusts=cutree(a$tree_row, 5)

kh29SubNames=kh29Aln$nam[ids]
names(kh29Clusts)=kh29SubNames

###################################
#
# Q: now plot these variants on the world map

#load ecotypes
#read in ecotype maps
ecotypes=read.csv('data/ecotype/ecotypes.csv', header=T, encoding="Latin-1")

#Do kh17 first

group1=sort(as.numeric(sapply(strsplit(kh17SubNames[which(kh17Clusts==1)], "\\|"), function(i){i[4]})))
group1Eco=ecotypes[which(ecotypes$pk %in% group1),]

group2=sort(as.numeric(sapply(strsplit(kh17SubNames[which(kh17Clusts==2)], "\\|"), function(i){i[4]})))
group2Eco=ecotypes[which(ecotypes$pk %in% group2),]

group3=sort(as.numeric(sapply(strsplit(kh17SubNames[which(kh17Clusts==3)], "\\|"), function(i){i[4]})))
group3Eco=ecotypes[which(ecotypes$pk %in% group3),]

group4=sort(as.numeric(sapply(strsplit(kh17SubNames[which(kh17Clusts==4)], "\\|"), function(i){i[4]})))
group4Eco=ecotypes[which(ecotypes$pk %in% group4),]

group5=sort(as.numeric(sapply(strsplit(kh17SubNames[which(kh17Clusts==5)], "\\|"), function(i){i[4]})))
group5Eco=ecotypes[which(ecotypes$pk %in% group5),]

latitude=list(group1Eco$latitude, group2Eco$latitude,group3Eco$latitude,group4Eco$latitude,group5Eco$latitude)
dev.off()
boxplot(latitude) #group 5 is the big conserved group; group 2 is the other group; group 1 is between


mergeEco=rbind(group1Eco, group2Eco, group3Eco, group4Eco, group5Eco)
mergeEco$groupolour=c(rep("lightblue", length(group1)),
                      rep("blue", length(group2)),
                      rep("grey", length(group3)),
                      rep("pink", length(group4)),
                      rep("red", length(group5)))

plot(mergeEco$latitude, mergeEco$longitude, col=mergeEco$groupolour, pch=20)

map(database="world", xlim=c(min(mergeEco$longitude, na.rm=T), max(mergeEco$longitude, na.rm=T)), 
    ylim=c(min(mergeEco$latitude, na.rm=T), max(mergeEco$latitude, na.rm=T))) 

# marking points on map 
points(x = mergeEco$longitude, y = mergeEco$latitude, col = mergeEco$groupolour, pch=19, cex=0.5)


##zoom in on Europe
map(database="world", xlim=c(-20, 100), 
    ylim=c(30, 65) )

# marking points on map 
points(x = mergeEco$longitude, y = mergeEco$latitude, col = mergeEco$groupolour, pch=19, cex=0.5)


#####################
# How does this correlate with ecotype flowering time?
ecotypeFT=read.csv('data/ecotype/ecotypeFloweringTimes.csv', header=T)

lowerNames=sapply(paste(ecotypes$name), function(i){
  print(i)
  try(tolower(paste(i)))})

nameMap=cbind(ecotypes$name, lowerNames)

easyIDs=lapply(ecotypeFT$Accession, function(i){
  a=tolower(i)
  a=grep(a, lowerNames)
  if(length(a)==0){print(i)}else{a}
})

#add appropriate name
accessionMatched=sapply(easyIDs, function(i){lowerNames[i[1]]})

ecotypeFT$accessionLower=accessionMatched


#which of these have good kh17 sequences?
lowerNamesMerge=sapply(paste(mergeEco$name), function(i){
  print(i)
  try(tolower(paste(i)))})

col=apply(ecotypeFT, 1, function(i){
  temp=i['accessionLower']
  me_id=which(lowerNamesMerge==temp)
  if(length(me_id>0)){
    mergeEco$groupolour[me_id]
  }else{'white'}
})


#############################
# What is the association between flowering time, flowering stdev and kh17?

ecotypeFT$latitude=sapply(ecotypeFT$accessionLower, function(i){
  mergeEco$latitude[which(lowerNamesMerge==i)]
})

ecotypeFT$longitude=sapply(ecotypeFT$accessionLower, function(i){
  mergeEco$longitude[which(lowerNamesMerge==i)]
})

plot(ecotypeFT$mean.ft, ecotypeFT$se.ft, col=col)
incl=which(col!='white')
ecokh17=ecotypeFT[incl,]
rownames(ecokh17)=ecokh17$Accession
orderFT=order(as.numeric(ecokh17$longitude))

pdf('plots/Ecotypes/Fig_NaturalEcotypeFloweringTime.pdf', width=8, height=14)
par(mfrow=c(3, 1))

##zoom in on Europe
map(database="world", xlim=c(-20, 100), 
    ylim=c(30, 65) )

# marking points on map 
points(x = ecokh17$longitude, y = ecokh17$latitude, labels=ecokh17$Accessions, col = col[incl], pch=19, cex=1)
points(x = mergeEco$longitude, y = mergeEco$latitude, col = mergeEco$groupolour, cex=0.5)

barplot(ecokh17$mean.ft[orderFT], col=col[incl[orderFT]], names=ecokh17$Accession[orderFT], las=2, ylab='flowering time mean (d)')
barplot(ecokh17$se.ft[orderFT], col=col[incl[orderFT]], names=ecokh17$Accession[orderFT], las=2, ylab='flowering time standard error (d)')


dev.off()


##########################################################
##make map-related plots for KH29

group1=sort(as.numeric(sapply(strsplit(kh29SubNames[which(kh29Clusts==1)], "\\|"), function(i){i[4]})))
group1Eco=ecotypes[which(ecotypes$pk %in% group1),]

group2=sort(as.numeric(sapply(strsplit(kh29SubNames[which(kh29Clusts==2)], "\\|"), function(i){i[4]})))
group2Eco=ecotypes[which(ecotypes$pk %in% group2),]

group3=sort(as.numeric(sapply(strsplit(kh29SubNames[which(kh29Clusts==3)], "\\|"), function(i){i[4]})))
group3Eco=ecotypes[which(ecotypes$pk %in% group3),]

group4=sort(as.numeric(sapply(strsplit(kh29SubNames[which(kh29Clusts==4)], "\\|"), function(i){i[4]})))
group4Eco=ecotypes[which(ecotypes$pk %in% group4),]

group5=sort(as.numeric(sapply(strsplit(kh29SubNames[which(kh29Clusts==5)], "\\|"), function(i){i[4]})))
group5Eco=ecotypes[which(ecotypes$pk %in% group5),]

latitude=list(group1Eco$latitude, group2Eco$latitude,group3Eco$latitude,group4Eco$latitude,group5Eco$latitude)
dev.off()
boxplot(latitude) #group 5 is the big conserved group; group 2 is the other group; group 1 is between


mergeEco=rbind(group1Eco, group2Eco, group3Eco, group4Eco, group5Eco)
mergeEco$groupolour=c(rep("blue", length(group1)),
                      rep("grey", length(group2)),
                      rep("grey", length(group3)),
                      rep("pink", length(group4)),
                      rep("red", length(group5)))

plot(mergeEco$latitude, mergeEco$longitude, col=mergeEco$groupolour, pch=20)

map(database="world", xlim=c(min(mergeEco$longitude, na.rm=T), max(mergeEco$longitude, na.rm=T)), 
    ylim=c(min(mergeEco$latitude, na.rm=T), max(mergeEco$latitude, na.rm=T))) 

# marking points on map 
points(x = mergeEco$longitude, y = mergeEco$latitude, col = mergeEco$groupolour, pch=19, cex=0.5)


##zoom in on Europe
map(database="world", xlim=c(-20, 100), 
    ylim=c(30, 65) )

# marking points on map 
points(x = mergeEco$longitude, y = mergeEco$latitude, col = mergeEco$groupolour, pch=19, cex=0.5)


#####################
# How does this correlate with ecotype flowering time?
ecotypeFT=read.csv('data/ecotype/ecotypeFloweringTimes.csv', header=T)

lowerNames=sapply(paste(ecotypes$name), function(i){
  print(i)
  try(tolower(paste(i)))})

easyIDs=lapply(ecotypeFT$Accession, function(i){
  a=tolower(i)
  a=grep(a, lowerNames)
  if(length(a)==0){print(i)}else{a}
})

#add appropriate name
accessionMatched=sapply(easyIDs, function(i){lowerNames[i[1]]})

ecotypeFT$accessionLower=accessionMatched

#which of these have good kh29 sequences?
lowerNamesMerge=sapply(paste(mergeEco$name), function(i){
  print(i)
  try(tolower(paste(i)))})

col=apply(ecotypeFT, 1, function(i){
  temp=i['accessionLower']
  me_id=which(lowerNamesMerge==temp)
  if(length(me_id>0)){
    mergeEco$groupolour[me_id]
  }else{'white'}
})



#############################
# What is the association between flowering time, flowering stdev and kh29?

ecotypeFT$latitude=sapply(ecotypeFT$accessionLower, function(i){
  mergeEco$latitude[grep(i, lowerNamesMerge)]
})

ecotypeFT$longitude=sapply(ecotypeFT$accessionLower, function(i){
  mergeEco$longitude[which(lowerNamesMerge==i)]
})

plot(ecotypeFT$mean.ft, ecotypeFT$se.ft, col=col)
incl=which(col!='white')
ecokh29=ecotypeFT[incl,]
rownames(ecokh29)=ecokh29$Accession
orderFT=order(as.numeric(ecokh29$longitude))

pdf('plots/Ecotypes/Fig_NaturalEcotypeFloweringTime_kh29.pdf', width=8, height=14)
par(mfrow=c(3, 1))

##zoom in on Europe
map(database="world", xlim=c(-20, 100), 
    ylim=c(30, 65) )

# marking points on map 
points(x = ecokh29$longitude, y = ecokh29$latitude, labels=ecokh29$Accessions, col = col[incl], pch=19, cex=1)
points(x = mergeEco$longitude, y = mergeEco$latitude, col = mergeEco$groupolour, cex=0.5)

barplot(ecokh29$mean.ft[orderFT], col=col[incl[orderFT]], names=ecokh29$Accession[orderFT], las=2, ylab='flowering time mean (d)')
barplot(ecokh29$se.ft[orderFT], col=col[incl[orderFT]], names=ecokh29$Accession[orderFT], las=2, ylab='flowering time standard error (d)')


dev.off()


