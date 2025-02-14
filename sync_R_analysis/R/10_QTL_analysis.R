#QTL mapping using the map by sequencing =================

#read in map with additional traits 
trait_table <- read.cross("csv", "./", "FULL_table_of_traits_map.csv", estimate.map=FALSE, crosstype = "riself", genotypes=c("W","X","H","U"))

summary(trait_table)
plotMissing(trait_table)

#want to look at individuals and markers with lots of missing information 
plot(ntyped(trait_table, "ind"), ylab="No.markers", xlab = "No. individuals", main="No. genotypes by individual")
plot(ntyped(trait_table, "mar"), ylab="No. individuals",xlab = "No. markers", main="No. genotypes by marker")

# In this case with such high numbers of markers its still not missing in more than 50% so will remove just one individual for this 
mapthis <- subset(trait_table, ind=(ntyped(trait_table)>600))
summary(mapthis)

#with the markers there are some markers that are missing in lots of individual and due to the amount of markers it is good to remove 

#To omit the markers with lots of missing data, we first need to identify the names of the markers, then use drop.markers().
nt.bymar <- ntyped(mapthis, "mar")
todrop <- names(nt.bymar[nt.bymar < 40])
mapthis <- drop.markers(mapthis, todrop)

summary(mapthis)

#drop douplicate markers
dup <- findDupMarkers(mapthis, exact.only=TRUE, adjacent.only = TRUE)

mapthis <- drop.markers(mapthis, unlist(dup))
summary(mapthis)

# filter genotypes not seggregating at 1:1
gt <- geno.table(mapthis)
summary(mapthis)

gt[gt$P.value < 0.05/totmar(mapthis),]
todrop <- rownames(gt[gt$P.value < 1e-10,])

mapthis <- drop.markers(mapthis, todrop)
summary(mapthis)

#to look at the genitype frequency 
g <- pull.geno(mapthis)

gfreq <- apply(g, 1, function(a) table(factor(a, levels=1:2)))
gfreq <- t(t(gfreq) / colSums(gfreq))

# filter duplicated records

plot(gfreq[1,], ylab="Genotype frequency", main=c("Ws", "Tnz")[1], ylim=c(0,1))
plot(gfreq[2,], ylab="Genotype frequency", main=c("Ws", "Tnz")[2], ylim=c(0,1))

# recombination fraction
mapthis <- est.rf(mapthis)
checkAlleles(mapthis, threshold=5)
rf <- pull.rf(mapthis)
lod <- pull.rf(mapthis, what="lod")
plotRF(mapthis)
# plot recombination fration
plot(as.numeric(rf), as.numeric(lod), xlab="Recombination fraction", ylab="LOD score")

#plot missingness again 
plotMissing(mapthis, main="missing data (after filtering)")

# look at the gentic map =============

#estimate the genetic distance between markers 
map <- est.map(mapthis, error.prob=0.0001, map.function="kosambi")

mapthis <- replace.map(mapthis, map)

plot.map(mapthis)

summary(mapthis)
#save new map .... 

write.cross(mapthis, format=c("csv"), filestem="filtered_map")


#Using map perform qtl anaylsis ===============

set.seed(1234)
sug <- calc.genoprob(mapthis, step=0, map.function="kosambi", error.prob = 0.0001)

out_hk <- scanone(sug, pheno.col = 2:990, method = "hk")

# only run this thresh_HK onces as time consuming for permutation
thresh_HK <- scanone(sug, pheno.col = 2:990, n.perm = 1000, method = "hk")
threshold_5_HK <- summary(thresh_HK, alpha = 0.05)

#save Rdata so dont have to keep running through
save.image("QTL_maping_lod.RData")

#save the outhk with lod values for each trait 
write.csv(out_hk, "QTL_output_LODscores.csv")
#save the permutation threshold as well 
write.csv(threshold_5_HK, "threshold_QTLs.csv")

#    Fliter even more for things that interact with developmental traits ================


#look at origional scores but only look at things which interact with the developmental traits make 0 and 1 for if there is a significant QTL present

# Loop through each column and if the lod score in out_hk is great than the threshold then it means there is a QTL present so replaces with a 1 and if not then a 0
for (i in 3:ncol(out_hk)) {
  out_hk[,i] <- ifelse(out_hk[,i] > threshold_5_HK[(i-2)], 1, 0)
}

triat <- out_hk[,-c(1,2)]
triat<- triat[,colSums(triat, na.rm = TRUE) != 0]
#and remove rows that contain zero too
triat <- triat[rowSums(triat, na.rm = TRUE) != 0 ,]

#transpose
triat <- t(triat)

#if the colsums for the rows containing developmental traits is > 0 keep (development taits are in rows: )   
filtered_develop_traits <- triat[,colSums(triat[273:278,], na.rm = TRUE) != 0]

#Then want to remove any rows with the rowsum 0 
filtered_develop_traits <- filtered_develop_traits[rowSums(filtered_develop_traits, na.rm = TRUE) != 0 ,]

#save the table of triats and regions of QTLs
write.csv(filtered_develop_traits, "QTLS_developmental_traits.csv")


network <- graph_from_biadjacency_matrix(filtered_develop_traits, weighted = NULL, directed = FALSE, mode = "out")

#what to colour by SNP (square), and functional trait (circle blue ) and developmental trait (circile darkblue)

#need vector that is 1111, 222, 3333333, 

# the node are set, this will seperate the functional traits (), developmental traits(5) and all the snps 
catogries <- c(rep(1,66),rep(2,6), rep(3,86))

# define color and shape mappings.
col <- c("steelblue","darkblue", "orange")
shape <- c("circle","circle", "square")
size <- c(15,15,3)
label_size <- c(0.8, 0.8, 0.8)



plot(network,
     vertex.color = col[catogries],
     vertex.shape = shape[catogries], lwd = 2, 
     vertex.size= size[catogries],  
     vertex.label.cex=label_size[catogries],
     vertex.label.dist = 1,
     vertex.label.color = "black", margin = c(0,0,0,0), rescale = T, )

tkplot(network,vertex.color = col[catogries],
       vertex.shape = shape[catogries], lwd = 2, 
       vertex.size= size[catogries],  
       vertex.label.cex=label_size[catogries],
       vertex.label.dist = 1,
       vertex.label.color = "black", margin = c(0,0,0,0), rescale = T)


#save.image("C:/Users/sl1407/OneDrive - University of York/Desktop/POSTDOC/TOPCOUNT/RIL_WxT_SD_LD/Final_scripts/QTL_maping_developmental.RData")




# get pvalue tables=========

# do pvalue at 10% 
threshold_10_HK <-  summary(thresh_HK, alpha = 0.10)


#reload outhk with lod scores
set.seed(1234)
out_hk <- scanone(sug, pheno.col = 2:990, method = "hk")

summary_lod_bypheno_pvalues <- summary(out_hk, perms = thresh_HK, format="allpheno", pvalues = TRUE)

#remove all the columns that contrain p.vales 
oddcols <- seq(3,1980, by = 2)

summary_lod_bypheno <- summary_lod_bypheno_pvalues[,oddcols]
summary_lod_bypheno <- cbind(summary_lod_bypheno_pvalues[,1:2],summary_lod_bypheno)

# Loop through each column and if the lod score in out_hk is great than the threshold then it means there is a QTL present so replaces with a 1 and if not then a 0

for (i in 3:ncol(summary_lod_bypheno)) {
  summary_lod_bypheno[,i] <- ifelse(summary_lod_bypheno[,i] > threshold_10_HK[(i-2)], 1, 0)
}

# dataframe of each of the pvalues with the same headings as the dataframe with the QTL presence and absence 
evencols <- seq(4,1980, by = 2)

summary_lod_bypheno_pvalues <- summary_lod_bypheno_pvalues[,evencols]

names_oftraits <- colnames(summary_lod_bypheno[,-c(1,2)])

colnames(summary_lod_bypheno_pvalues) <- names_oftraits


#replacing the 1s with the pvalues 
summary_lod_bypheno <- summary_lod_bypheno[,-c(1,2)]

#first must replace all the 0s in summary_lod_bypheno with dummy number (200)
summary_lod_bypheno[summary_lod_bypheno == 0] <- 200

summary_lod_bypheno[summary_lod_bypheno == 1] <- summary_lod_bypheno_pvalues[summary_lod_bypheno == 1]

#replace dummy number with NA
summary_lod_bypheno[summary_lod_bypheno == 200] <- NA


#remove columns and rows that total 0 i.e no QTLs 
summary_lod_bypheno<- summary_lod_bypheno[,colSums(summary_lod_bypheno, na.rm = TRUE) != 0]
#and remove rows (SNP markers) that are zero too i.e no QTL present 
summary_lod_bypheno <- summary_lod_bypheno[rowSums(summary_lod_bypheno, na.rm = TRUE) != 0 ,]

#save files 

write.csv(summary_lod_bypheno, "DEG_ws_tnz_with_pavlue_Thresh10.csv")




