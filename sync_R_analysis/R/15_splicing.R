###########
# The aim of this section is to explore alternative splicing in FLC family members
# Specifically, we will look to see if splice variants are linked to KH alleles
##########


#######
# For the important lines, look at the splice variants over gene expression


flcfamily=c("AT5G10140", "AT1G77080", "AT5G65050", "AT5G65060", "AT5G65070", "AT5G65080")
chr=      c(5, 1, 5, 5, 5, 5)
starts=   c(3173382, 28955522, 25982254, 25987309, 25992243, 25997594)
stops=    c(3179448, 28960111, 25986435, 25991412, 25996174, 26002330)

getSplicesInRange <- function(chr, start, stop, sj){
  temp=sj[which(sj[,1]==paste("Chr", chr, sep="")),]
  temp=temp[which(temp[,2]>=start),]
  temp[which(temp[,3]<=stop),]
}


flcGene=1
line=1
importantLines_T=unique(gsub("_C", "_T", list.dirs(path='data\\BAMs\\mapping_output\\mapping_output', full.names = F, recursive=F)))
importantLines_C=unique(gsub("_T", "_C", list.dirs(path='data\\BAMs\\mapping_output\\mapping_output', full.names = F, recursive=F)))

splicesInKeyLines=lapply(1:length(importantLines_C), function(line){
  
  #T
  fileT=read.table(paste('data\\BAMs\\mapping_output\\mapping_output\\', importantLines_T[line], "\\SJ.out.tab", sep=""))
  
  #C
  fileC=read.table(paste('data\\BAMs\\mapping_output\\mapping_output\\', importantLines_C[line], "\\SJ.out.tab", sep=""))
  
  lapply(1:length(flcfamily), function(flcGene){
    
    fileT_sj=getSplicesInRange(chr[flcGene], starts[flcGene], stops[flcGene], fileT)
    fileC_sj=getSplicesInRange(chr[flcGene], starts[flcGene], stops[flcGene], fileC)
    
    list(T_sj=fileT_sj, C_sj=fileC_sj)
    
  })
})

#find unique splice junctions
uniqueSJ=lapply(1:length(flcfamily), function(flcGene){
  unique(unlist(lapply(1:length(importantLines_C), function(line){
    flc=splicesInKeyLines[[line]][[flcGene]]
    a=paste(flc[[1]][,1], flc[[1]][,2], flc[[1]][,3])
    b=paste(flc[[2]][,1], flc[[2]][,2], flc[[2]][,3])
    unique(c(a,b))
  })))
})

uniqueTablesMergedPerGene=lapply(1:length(uniqueSJ), function(flcGene){
  line_sj=lapply(1:length(importantLines_C), function(line){
    print("start")
    print(paste(line, flcGene))
    
    
    flc_t=splicesInKeyLines[[line]][[flcGene]][[1]]
    flc_c=splicesInKeyLines[[line]][[flcGene]][[2]]
    
    
    sj_id_t=paste(flc_t[,1], flc_t[,2], flc_t[,3])
    sj_id_c=paste(flc_c[,1], flc_c[,2], flc_c[,3])
    
    
    tC=sapply(uniqueSJ[[flcGene]], function(sj){
      if(sj %in% sj_id_t){
        flc_t[which(sj_id_t==sj),7][1]
        
      }else{0}
    })
    
    cC=sapply(uniqueSJ[[flcGene]], function(sj){
      if(sj %in% sj_id_c){
        flc_c[which(sj_id_c==sj),7][1]
        
      }else{0}
    })
    
    cbind(tC, cC)
  })
  
  print(dim(line_sj[[1]]))
  do.call(cbind, line_sj)
  
})

#########################################
# Now, make a nice heatmap of splice junctions, split by variants in KH genes

#normalise separately for each gene
uniqueTablesMergedPerGeneNorm=lapply(uniqueTablesMergedPerGene, function(tab){
  apply(tab, 2, function(j){j/sum(j)})
})

#make a nice splice map
map=do.call(rbind, uniqueTablesMergedPerGeneNorm)
map_count=do.call(rbind, uniqueTablesMergedPerGene)
#column names
colnames(map)=c(rbind(importantLines_C,importantLines_T))
colnames(map_count)=colnames(map)
annotCol=data.frame(row.names = colnames(map), 
                    tc=rep(c("C", "T"), length(importantLines_C)))

write.csv(map_count, 'data/BAMs/spliceJunctionsWithinGene.csv')

annotCol$kh17="W"
annotCol$kh29='W'


########################
########################################
# Look at the variants

#######
# let us look at the actual variants present

vcf=read.table(file = "data/BAMs/all_Control_combined.vcf", header=F)

colVCF=read.table(file = "data/BAMs/all_Control_combined.vcf", header=F, comment.char="", skip = 44, nrows=1)

colnames(vcf)=colVCF

#only include kh17 variants
chr=3
start=13490828
stop=13493834


vcf=vcf[which(vcf[,1]==paste("Chr", chr, sep='')),]

vcf=vcf[which(as.numeric(vcf[,2])>=start),]

vcf=vcf[which(as.numeric(vcf[,2])<=stop),]

cleaned=apply(vcf, 1, function(i){
  idsW=grep("0/0", i)
  idsH=c(grep("0/1", i),grep("1/0", i))
  idsT=grep("1/1", i)
  idsM=grep("./.", i)
  temp=i
  temp[idsM]="M"
  temp[idsW]="W"
  temp[idsH]="H"
  temp[idsT]="T"
  temp
})

cleaned=t(cleaned)


#merge variants for split lanes
l3=grep("L3", colnames(cleaned))
l4=grep("L4", colnames(cleaned))

replace=sapply(1:length(l3), function(i){
  a=cleaned[,l3[i]]
  b=cleaned[,l4[i]]
  
  temp=a
  temp[which(a=="M")]=b[which(a=="M")]
  temp[which(a=="H" | b=="H")]="H"
  temp
})

cleaned[,l3]=replace
cleaned=cleaned[,-l4]

variantDist=apply(cleaned[,-c(1:9)], 1, function(i){table(factor(i, levels=c("W", "T", "H", "M")))})

#remove variants where either T or W are zero
cleaned=cleaned[which(variantDist[1,]!=0 & variantDist[2,]!=0),]

variantDist=apply(cleaned[,-c(1:9)], 1, function(i){table(factor(i, levels=c("W", "T", "H", "M")))})

#replace T and W so that W is the most common variant

cleaned=t(apply(cleaned, 1, function(i){
  if(length(which(i=="T"))>length(which(i=="W"))){
    
    temp=i
    temp[which(i=='W')]='T'
    temp[which(i=='T')]='W'
    temp
  }else{i}
}))


variantDist=apply(cleaned[,-c(1:9)], 1, function(i){table(factor(i, levels=c("W", "T", "H", "M")))})

annotCols=data.frame(t(cleaned[,-c(1:9)]))
rownames(annotCols)=gsub("_L3_results", "", rownames(annotCols))

countType=apply(annotCols, 1, function(i){
  table(factor(i, levels=c('W', 'T', 'H', 'M')))
})

classification_kh17=apply(countType, 2, function(i){
  if(i['H']>(1/3)*sum(i[1:3])){'H'}else{
    
    c('W', 'T', 'H')[which.max(i[1:3])]
  }
})





vcf=read.table(file = "data/BAMs/all_Control_combined.vcf", header=F)

colVCF=read.table(file = "data/BAMs/all_Control_combined.vcf", header=F, comment.char="", skip = 44, nrows=1)

colnames(vcf)=colVCF

#only include kh29 variants
chr=5
#22725238 - 22728473
start=22725238
stop=22728473

vcf=vcf[which(vcf[,1]==paste("Chr", chr, sep='')),]

vcf=vcf[which(as.numeric(vcf[,2])>=start),]

vcf=vcf[which(as.numeric(vcf[,2])<=stop),]

cleaned=apply(vcf, 1, function(i){
  idsW=grep("0/0", i)
  idsH=c(grep("0/1", i),grep("1/0", i))
  idsT=grep("1/1", i)
  idsM=grep("./.", i)
  temp=i
  temp[idsM]="M"
  temp[idsW]="W"
  temp[idsH]="H"
  temp[idsT]="T"
  temp
})

cleaned=t(cleaned)


#merge variants for split lanes
l3=grep("L3", colnames(cleaned))
l4=grep("L4", colnames(cleaned))

replace=sapply(1:length(l3), function(i){
  a=cleaned[,l3[i]]
  b=cleaned[,l4[i]]
  
  temp=a
  temp[which(a=="M")]=b[which(a=="M")]
  temp[which(a=="H" | b=="H")]="H"
  temp
})

cleaned[,l3]=replace
cleaned=cleaned[,-l4]

variantDist=apply(cleaned[,-c(1:9)], 1, function(i){table(factor(i, levels=c("W", "T", "H", "M")))})

#remove variants where either T or W are zero
cleaned=cleaned[which(variantDist[1,]!=0 & variantDist[2,]!=0),]

variantDist=apply(cleaned[,-c(1:9)], 1, function(i){table(factor(i, levels=c("W", "T", "H", "M")))})

#replace T and W so that W is the most common variant

cleaned=t(apply(cleaned, 1, function(i){
  if(length(which(i=="T"))>length(which(i=="W"))){
    
    temp=i
    temp[which(i=='W')]='T'
    temp[which(i=='T')]='W'
    temp
  }else{i}
}))


variantDist=apply(cleaned[,-c(1:9)], 1, function(i){table(factor(i, levels=c("W", "T", "H", "M")))})

annotCols=data.frame(t(cleaned[,-c(1:9)]))
rownames(annotCols)=gsub("_L3_results", "", rownames(annotCols))

countType=apply(annotCols, 1, function(i){
  table(factor(i, levels=c('W', 'T', 'H', 'M')))
})


classification_kh29=apply(countType, 2, function(i){
  if(i['H']>(1/3)*sum(i[1:3])){'H'}else{
    
    c('W', 'T', 'H')[which.max(i[1:3])]
  }
})

###let us add the _T to the two classification variables
class_kh17=c(classification_kh17, classification_kh17[-c(1:4)])
names(class_kh17)=c(names(classification_kh17), 
                    gsub('_C', '_T', names(classification_kh17[-c(1:4)])))

class_kh29=c(classification_kh29, classification_kh29[-c(1:4)])
names(class_kh29)=c(names(classification_kh29), 
                    gsub('_C', '_T', names(classification_kh29[-c(1:4)])))



########################


annotCol[names(class_kh17)[which(class_kh17=='T')],'kh17']="T"
annotCol[names(class_kh17)[which(class_kh17=='H')],'kh17']="H"

annotCol[names(class_kh29)[which(class_kh29=='T')],'kh29']="T"
annotCol[names(class_kh29)[which(class_kh29=='H')],'kh29']="H"

#merge rare rows into a rareRowCol
thresh=0.001*length(rowSums(map))
idsCommon=which(rowSums(map)>thresh)
sumUncommon=colSums(map[-idsCommon,])/length(which(rowSums(map)<=thresh))
map=rbind(map[idsCommon,], sumUncommon)

#how much do the sequences correlate with each other?
crosscor=apply(annotCols, 1, function(i){
  apply(annotCols, 1, function(j){
    length(which(i==j))
  })
})

a=pheatmap(crosscor)
kh17Clusts=cutree(a$tree_col, 3)
kh17Clusts2=kh17Clusts
names(kh17Clusts2)=gsub("_C", "_T", names(kh17Clusts))
kh17Clusts=c(kh17Clusts, kh17Clusts2)
annotCol$variants='W'
annotCol[unique(names(kh17Clusts)[which(kh17Clusts==1)]), "variants"]='T'
annotCol[unique(names(kh17Clusts)[which(kh17Clusts==3)]), "variants"]='H'
annotCol$variants=factor(annotCol$variants)
print('about to do pheatmap')
#if(g==5 | g==6){
#  print(dim(map))
#  rm=which(apply(map, 2, function(i){sd(i)})==0)
#  print(rm)
#  map=map[,-rm]
#}

annot_row=data.frame(row.names=rownames(map))
annot_row$gene=as.factor(sapply(rownames(map), function(i){
  if(i=='sumUncommon'){
    'sum'
  }else{
    flcfamily[which(sapply(uniqueTablesMergedPerGene, function(j){
      i %in% rownames(j)
    }))]
  }
}))


annot_row$gene=as.factor(annot_row$gene)

annot_colours=list(
  kh29=c("H"='violet', "T"='blue', "W"="lightgrey"),
  kh17=c("H"='violet', "T"='blue', "W"="lightgrey"),
  tc=c("C"="black", "T"="lightyellow"),
  gene=c("AT1G77080"="lightgreen", "AT5G10140"="darkgreen", "AT5G65050"="lightblue", "AT5G65060"="cornflowerblue", "sum"="white")
)

variantSet=pheatmap(map, scale='none', annotation_colors=annot_colours, annotation_row=data.frame(annot_row), annotation_col = annotCol, main="combined", cutree_cols=2)
cl=cutree(variantSet$tree_col,2)
#grab most different ones and make condensed plot
id_variants=which(sapply(1:dim(map)[1], function(i){
  a=map[i, names(cl)[which(cl==1)]]
  b=map[i, names(cl)[which(cl==2)]]
  t.test(a, b)$p.value<0.05
}))
variantSet=pheatmap(map[id_variants,], scale='row',labels_col=NA, color = colorRampPalette(rev(brewer.pal(n = 11, name =
                                                                                                            "RdBu")))(100), annotation_colors=annot_colours, annotation_row=data.frame(annot_row), annotation_col = annotCol, main="combined", cutree_cols=2)

annotCol$KH=paste(annotCol$kh17, annotCol$kh29)

annot_colours=list(
  KH=c("T T"="black", "T H"="mediumblue", "T W"="lightblue", 'H T'='mediumblue', 'W T'='lightblue', 'H W'='pink', 'W H'='pink', 'W W'='white'),
  
  # KH=c("T T"="black", "T H"="darkblue", "T W"="cornflowerblue", 'H T'='darkmagenta', 'W T'='violet', 'H W'='lightblue', 'W H'='pink', 'W W'='white'),
  kh29=c("H"='violet', "T"='blue', "W"="lightgrey"),
  kh17=c("H"='violet', "T"='blue', "W"="lightgrey"),
  tc=c("C"="black", "T"="lightyellow"),
  gene=c("AT1G77080"="yellow", "AT5G65050"="darkgreen", "AT5G65060"="lightgreen", "sum"="white")
)



variantSet=pheatmap(map[id_variants,], scale='row',labels_col='line', color = colorRampPalette(rev(brewer.pal(n = 11, name =
                                                                                                                "RdBu")))(100), annotation_colors=annot_colours, annotation_row=data.frame(annot_row), annotation_col = annotCol[,c('tc','KH')], main="combined", cutree_cols=2, file='plots/Splice/Fig_KH_spliceHeatmap.pdf')

dev.off()


#############
# Now analyse chimeric transcripts
###############################################
###########Explore super weird idea.  Are there MAF2/MAF3 merged transcripts?



flcfamilyImp=c("AT5G65050", "AT5G65060")
chrImp=      c(5, 5)
startsImp=   c(25982254, 25987309)
stopsImp=    c(25986435, 25991412)

getSplicesInRangeMergedGenes <- function(chr, start1, start2, stop1, stop2, sj){
  temp=sj[which(sj[,1]==paste("Chr", chr, sep="")),]
  temp[which(temp[,2]>=start1 & temp[,2]<=stop1 & temp[,3]>=start2 & temp[,3]<=stop2),]
}


flcGene=1
line=1
#importantLines_T=unique(gsub("_C", "_T", list.dirs(path='.', full.names = F, recursive=F)))
#importantLines_C=unique(gsub("_T", "_C", list.dirs(path='.', full.names = F, recursive=F)))

splicesInKeyLinesMerged=lapply(1:length(importantLines_C), function(line){
  
  #T
  fileT=read.table(paste('data\\BAMs\\mapping_output\\mapping_output\\', importantLines_T[line], "\\SJ.out.tab", sep=""))
  
  #C
  fileC=read.table(paste('data\\BAMs\\mapping_output\\mapping_output\\', importantLines_C[line], "\\SJ.out.tab", sep=""))
  ##key change
  
  fileT_sj=getSplicesInRangeMergedGenes(chrImp[1], startsImp[1], startsImp[2], stopsImp[1], stopsImp[2], fileT)
  fileC_sj=getSplicesInRangeMergedGenes(chrImp[1], startsImp[1], startsImp[2], stopsImp[1], stopsImp[2], fileC)
  
  list(T_sj=fileT_sj, C_sj=fileC_sj)
  
  
})


weirdSplices=lapply(splicesInKeyLinesMerged, function(i){
  temp1=data.frame(names=paste(i$C_sj[,1], i$C_sj[,2], i$C_sj[,3]), count=i$C_sj[,7])
  if(dim(i$C_sj)[1]>0){
    temp1$condition='C'
  }
  
  temp2=data.frame(names=paste(i$T_sj[,1], i$T_sj[,2], i$T_sj[,3]), count=i$T_sj[,7])
  if(dim(i$T_sj)[1]>0){
    temp2$condition='T'
  }
  
  print(dim(temp1))
  print(dim(temp2))
  temp=rbind(temp1, temp2)
  temp[which(temp$count>0),]
})

names(weirdSplices)=importantLines_C

uniqueWeirdSplices=unique(unlist(sapply(weirdSplices, function(i){i$names})))

annotCol[names(class_kh17)[which(class_kh17=='T')],'kh17']="T"
annotCol[names(class_kh17)[which(class_kh17=='H')],'kh17']="H"

annotCol[names(class_kh29)[which(class_kh29=='T')],'kh29']="T"
annotCol[names(class_kh29)[which(class_kh29=='H')],'kh29']="H"


weirdSplicesC=sapply(uniqueWeirdSplices, function(i){
  sapply(weirdSplices, function(j){
    if(length(which(j$names==i & j$condition=='C'))>0){
      j[which(j$names==i & j$condition=='C'),'count']
    }else{0}
  })
})

weirdSplicesT=sapply(uniqueWeirdSplices, function(i){
  sapply(weirdSplices, function(j){
    if(length(which(j$names==i & j$condition=='T'))>0){
      j[which(j$names==i & j$condition=='T'),'count']
    }else{0}
  })
})

weirdSplicesT_rowNameRight=weirdSplicesT
rownames(weirdSplicesT_rowNameRight)=gsub("_C", "_T", rownames(weirdSplicesT_rowNameRight))

weirdSplicesCombined=rbind(weirdSplicesC, weirdSplicesT_rowNameRight)

write.csv(weirdSplicesCombined, 'data/BAMs/weirdSplices.csv')
plot(rowSums(weirdSplicesC), rowSums(weirdSplicesT), xlab='cross-transcript splice')
kh11Names=rownames(annotCol)[which(annotCol$kh11=='T')]
points(rowSums(weirdSplicesC)[which(rownames(weirdSplicesC) %in% kh11Names)], rowSums(weirdSplicesT)[which(rownames(weirdSplicesC) %in% kh11Names)], pch=19, col='blue')

plot(rowSums(weirdSplicesC), rowSums(weirdSplicesT), xlab='cross-transcript splice')
kh29Names=rownames(annotCol)[which(annotCol$kh29=='T')]
points(rowSums(weirdSplicesC)[which(rownames(weirdSplicesC) %in% kh29Names)], rowSums(weirdSplicesT)[which(rownames(weirdSplicesC) %in% kh29Names)], pch=19, col='blue')

png(file='plots/Splice/Fig_barplotChimeras.png', width=16, height=5, res=400, units='in')

barplot(sort(rowSums(weirdSplicesC)+rowSums(weirdSplicesT)), ylim=c(0, 1300), ylab='# MAF2/MAF3 chimeric splice junctions', xlab='RILs', axisnames=F)
width=1.2
sapply(1:length(weirdSplicesC[,1]), function(i){
  nm=names(sort(rowSums(weirdSplicesC)+rowSums(weirdSplicesT)))[i]
  col=annot_colours$KH[annotCol[nm,'KH']]
  rect((i-1)*width, 1050, i*width, 1300, col=col)
})
dev.off()

png(file='plots/Splice/Fig_barplotChimeras_K29.png', width=16, height=5, res=400, units='in')

cols_k29=sapply(1:length(weirdSplicesC[,1]), function(i){
  nm=names(sort(rowSums(weirdSplicesC)+rowSums(weirdSplicesT)))[i]
  annot_colours$kh29[annotCol[nm,'kh29']]
})
barplot(sort(rowSums(weirdSplicesC)+rowSums(weirdSplicesT)), ylim=c(0, 1300), ylab='# MAF2/3 chimeric SJs', xlab='RILs', axisnames=F, col=cols_k29, cex.lab=1.6)
width=1.2
sapply(1:length(weirdSplicesC[,1]), function(i){
  nm=names(sort(rowSums(weirdSplicesC)+rowSums(weirdSplicesT)))[i]
  #  col1=annot_colours$kh11[annotCol[nm,'kh11']]
  col2=annot_colours$kh29[annotCol[nm,'kh29']]
  rect((i-1)*width, 1100, i*width, 1200, col=col2)
  #  rect((i-1)*width, 1200, i*width, 1300, col=col1)
})
legend(0, 1000, c('Tnz', 'Hetero', 'Ws-2'), title='KH29 genotype', fill=c('blue', 'violet', 'grey'))

dev.off()

