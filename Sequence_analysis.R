##### Analysis of High Temporal Resolution Land Use/Land Cover Trajectories
##### Supplementary material to article published in Land https://www.mdpi.com/journal/land
##### R script to carry out the analysis presented in the paper
##### JF Mas - Universidad Nacional Autónoma de México - jfmas@ciga.unam.mx



# Working directory and libraries
setwd("/home/jf/Dropbox/sequences/SupplementaryMaterials")
library(TraMineR)
library(raster)
library(fastcluster)


###################################################################################
## Read rasters of tome serie of land cover and covariates

## LUC maps for eleven dates
serie <- stack("serie.tif")
summary(serie)
# Distance from roads
dist <- raster("dist_roads.tif")
# Elevation
dem <- raster("dem.tif")
# Slope
slope <- raster("slope.tif")

## Plot the maps
plot(serie[[1]])
plot(dist)
plot(dem)
plot(slope)

####################################################################################
## Elaboration of a table which describes for each pixel the land category at each time step + value of covariates

## using rasterToPoints(r) to get coordinates
tab <- as.data.frame(rasterToPoints(stack(serie,dem,slope,dist)))
names(tab) <- c("x","y","c1986","c1989","c1992","c1995","c1998","c2001","c2004",
                "c2007","c2010","c2013","c2016","dem","slope","dist")
# A look at the first rows of the table:
head(tab)
# and a summary:
summary(tab)

## Create categorical data from covariates
tab$demK<-cut(tab$dem, c(400,700,1000,1200))
tab$slopeK<-cut(tab$slope, c(-1,2,6,102))
tab$distK<-cut(tab$dist, c(0,500,1000,2000,5000,10000))
## Two categories for distance and slope (binary variable)
tab$distK2<-cut(tab$dist, c(0,1000,10000))
tab$slopeK2<-cut(tab$slope, c(-1,6,102))

# Put an extra column with the concatened sequence
tab$sec <- as.vector(paste(tab[,3],tab[,4],tab[,5],tab[,6],tab[,7],
      tab[,8],tab[,9],tab[,10],tab[,11],tab[,12],tab[,13],sep="-"))

# Verify all the categories for all the time steps
clases <- c()
for (year in 3:13){
  clases <- unique(c(clases,tab[,year]))
  }
print(sort(clases))
# [1]  3  4 15 18 21 25


### Determine the number of different sequences
unique_traj <- sort(unique(tab$sec))
length(unique_traj)  # 1651 different trajectories

# Determine the frequency of these trajectories
tab_freq <- as.data.frame(table(tab$sec))

# Sort the sequences from more to less frequent
tab_freq2 <- tab_freq[order(tab_freq$Freq,decreasing = TRUE),] 
head(tab_freq2)

# Creation of a state sequence object to which names are assigned
# short state names (for printed output) and long state labels (for the graphics' legend). 
# Also define the alphabet (be sure that the order of the names and labels is the same of the alphabet.
# and a palette of colors for display

## State sequence object
alphabet <- c(3,4,15,18,21,25)
labels <- c("Forest", "Savanna","Pasture", "Agriculture","Mosaic Agri/Past","Others")                                   
short_labels <- c("F","Sav","Past","Agri","Mosaic","Others")
palette <- c("darkolivegreen","green","yellow","red","orange","grey")
tab.seq <- seqdef(tab, 3:13, alphabet = alphabet, states = short_labels,
                   cpal = palette, labels = labels)

#######################################################################################
## Visualize the sequence data set 

# Plot 10 sequences in the tab.seq sequence object (chosen to show diversity)
some_seq <- tab.seq[c(19,5,973,976,34,84,930,3893,993,995),]
seqiplot(some_seq, with.legend = T, border = T, main = "Some sequences")

# Plot all the sequences in the data set, sorted by states from start.
seqIplot(tab.seq, sortv = "from.start", with.legend = T, main = "Sequences 1985-2017")

# Plot the 10 most frequent sequences.
seqfplot(tab.seq, with.legend = T, main="Most common sequences")
 
#######################################################################################
## Explore the sequence data set by computing and visualizing descriptive statistics
 

# Compute and plot the state distributions by time step. 
# With border = NA, borders surrounding the bars are removed. 
seqdplot(tab.seq, with.legend = T, border = NA,main="Land cover (states) distribution", ylab="Proportion of study area")

# Compute and plot the transversal entropy index (Landscape entropy over time)
seqHtplot(tab.seq, main = "Entropy", ylab="Entropy index value",xlab=("Time"))

#Plot the sequence of modal states (dominant land cover) of the transversal state distributions.
seqmsplot(tab.seq, with.legend = T, main ="Most frequent land cover")

# Plot the mean time spent in each land cover category.
seqmtplot(tab.seq, with.legend = T, main = "Permanence", ylab="Number of 3 years periods")

### Longitudinal turbulence and entropy indices 
# Computed for each pixel over the time
tab$Turb <- seqST(tab.seq, norm=FALSE)
tab$Entrop <- seqient(tab.seq, norm=TRUE, base=exp(1))

# Generate rasters which represent these indices
xyt <- as.data.frame(cbind(tab$x,tab$y,tab$Turb))
names(xyt) <- c("x","y","t")
coordinates(xyt) <- ~ x + y
gridded(xyt) <- TRUE
raster_turb <- raster(xyt)
plot(raster_turb)

xye <- as.data.frame(cbind(tab$x,tab$y,tab$Entrop))
names(xye) <- c("x","y","e")
head(xye)
coordinates(xye) <- ~ x + y
gridded(xye) <- TRUE
raster_entrop <- raster(xye)
plot(raster_entrop)

# Calculate the correlation between both indices
cor(tab$Turb,tab$Entrop) # 0.9394874

## Computes the transition rates
tr_rates <- seqtrate(tab.seq)
print(tr_rates)

#######################################################################################
# Compute distances between sequences using different dissimilarity indices

## OM with substitution costs based on transition
## probabilities and indel set as half the maximum
## substitution cost
costs.tr <- seqcost(tab.seq, method = "TRATE",with.missing = FALSE)
print(costs.tr)
dist.om1 <- seqdist(tab.seq, method = "OM",indel = costs.tr$indel, sm = costs.tr$sm,with.missing = F)
dim(dist.om1)

### OM based on features
tab_state_features <- data.frame(state=c(10,5,3,3,3,1))
costs.gower <- seqcost(tab.seq, method = "FEATURES",with.missing = FALSE,state.features = tab_state_features)
print(costs.gower)

dist.om2 <- seqdist(tab.seq, method = "OM",indel = costs.gower$indel, sm = costs.gower$sm,with.missing = F)
dim(dist.om2)

## LCS
dist.lcs <- seqdist(tab.seq, method = "LCS")

## LCP
dist.lcp <- seqdist(tab.seq, method = "LCP") 

# Elaboration a typology of the trajectories: build a Ward hierarchical clustering
# of the sequences from the different distances and retrieve for each cell sequence the
# cluster membership of the 5 class solution. 

## Cluster based on OM transition rates
clusterward_om1 <- hclust(as.dist(dist.om1),method="ward.D")
plot(clusterward_om1)
cl_om1 <- cutree(clusterward_om1, k = 5)
tab$clusterom1 <- cl_om1
head(tab)

## Cluster based on OM features
clusterward_om2 <- hclust(as.dist(dist.om2),method="ward.D")
plot(clusterward_om2)
cl_om2 <- cutree(clusterward_om2, k = 5)
tab$clusterom2 <- cl_om2
head(tab)

## Cluster based on LCS
clusterward_lcs <- hclust(as.dist(dist.lcs),method="ward.D")
plot(clusterward_lcs)
cl_lcs <- cutree(clusterward_lcs, k = 5)
tab$clusterlcs <- cl_lcs
head(tab)

## Cluster based on LCP
clusterward_lcp <- hclust(as.dist(dist.lcp),method="ward.D")
plot(clusterward_lcp)
cl_lcp <- cutree(clusterward_lcp, k = 5)
tab$clusterlcp <- cl_lcp
head(tab)

# Plot all the sequences within each cluster para los 4 métodos
# OM1
seqIplot(tab.seq, group = tab$clusterom1, sortv = "from.start")
# OM2
seqIplot(tab.seq, group = tab$clusterom2, sortv = "from.start")
# LCS
seqIplot(tab.seq, group = tab$clusterlcs, sortv = "from.start")
# LCP
seqIplot(tab.seq, group = tab$clusterlcp, sortv = "from.start")


######### Plot clusters' spatial distribution

# Elaborate raster OM1
xyz <- as.data.frame(cbind(tab$x,tab$y,tab$clusterom1))
names(xyz) <- c("x","y","z")
xyz <- xyz[complete.cases(xyz), ]
coordinates(xyz) <- ~ x + y
gridded(xyz) <- TRUE
raster_com1 <- raster(xyz)
plot(raster_com1)

# Elaborate raster OM2
xyz <- as.data.frame(cbind(tab$x,tab$y,tab$clusterom2))
names(xyz) <- c("x","y","z")
xyz <- xyz[complete.cases(xyz), ]
coordinates(xyz) <- ~ x + y
gridded(xyz) <- TRUE
raster_com2 <- raster(xyz)
plot(raster_com2)

# Elaborate raster LCS
xyz <- as.data.frame(cbind(tab$x,tab$y,tab$clusterlcs))
names(xyz) <- c("x","y","z")
xyz <- xyz[complete.cases(xyz), ]
coordinates(xyz) <- ~ x + y
gridded(xyz) <- TRUE
raster_clcs <- raster(xyz)
plot(raster_clcs)

# Elaborate raster LCP
xyz <- as.data.frame(cbind(tab$x,tab$y,tab$clusterlcp))
names(xyz) <- c("x","y","z")
xyz <- xyz[complete.cases(xyz), ]
coordinates(xyz) <- ~ x + y
gridded(xyz) <- TRUE
raster_clcp <- raster(xyz)
plot(raster_clcp)

#######################################################################################
## Run discrepancy analyses to study how sequences are related to covariates

# Compute and test the share of discrepancy explained by different categories on covariates 
da1 <- dissassoc(dist.om1, group = tab$slopeK, R = 50)
print(da1$stat)
da2 <- dissassoc(dist.om1, group = tab$distK, R = 50)
print(da2$stat)
da3 <- dissassoc(dist.om1, group = tab$distK2, R = 50)
print(da3$stat)


# Selecting event subsequences:
# The analysis was restricted to sequences that exhibit the state Mosaic

tabe.seq <- seqecreate(tab.seq, use.labels = FALSE)
mosaic <- seqecontain(tabe.seq, event.list = c("Mosaic"))
mosaic_tab <- tab[mosaic,]
mosaic.seq <- tab.seq <- seqdef(mosaic_tab, 3:13, alphabet = alphabet, states = short_labels,
                                       cpal = palette, labels = labels)
mosaic.seqe <- seqecreate(mosaic.seq, use.labels = FALSE)

# Look for frequent event subsequences and plot the 10 most frequent ones.
fsubseq <- seqefsub(mosaic.seqe, pmin.support = 0.05)
head(fsubseq)
# 10 Most common subsequences
plot(fsubseq[1:10], col = "grey98")

# Determine the subsequences of transitions which best discriminate the groups as
# areas close and faraway from roads
discr1 <- seqecmpgroup(fsubseq, group = mosaic_tab$distK2)
plot(discr1[1:10],cex=1,cex.legend=1,legend.title="Distance",cex.lab=0.8, cex.axis = 0.8)
# areas with moderate vs steep slope
discr2 <- seqecmpgroup(fsubseq, group = mosaic_tab$slopeK2)
plot(discr2[1:10],cex=1,cex.legend=1,legend.title="Slope",cex.lab=0.8, cex.axis = 0.8)
# clusters of sequences
discr3 <- seqecmpgroup(fsubseq, group = mosaic_tab$clusterom1)
plot(discr3[1:10],cex=1,cex.legend=1,legend.title="Clusters OM1",cex.lab=0.8, cex.axis = 0.8)


