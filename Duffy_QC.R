setwd('/Users/jimblotter/Desktop/Grad_School/Data_Analysis/QC/') #remember to put the .str file there!
install.packages("genetics")
library("ape")
library("genetics")        #there is no package called ‘genetics’ = install.packages("genetics")
library("pegas")
library("seqinr")
library("ggplot2")
library("adegenet")

#Read in data
obj1 <- read.structure("eri_sor.str") #place cursor in console
# It will prompt for info:
#   genotypes = 246 (number of samples) GOT THIS FROM BARCODES FILE
#   markers = 1886 (number of loci)
#   column with labels for genotypes = 1 (the names of the samples)
#   column with population factor = 0 (if we were specifying populations a priori but we aren't)
#   other optional columns - just hit enter, there are none in our data
#   row with the marker names = 0 (We don't have a marker names row)
#   genotypes coded by a single row = Y (We have 2 lines per sample in the dataset? I THINK I JUST HAVE ONE)
indNames(obj1) # should return a list of the 267 sample names DOESNT
ploidy(obj1) # should return 2 since we gave it 2 alleles for each marker DOES

#Neighbor joining euclidian distance tree
D <- dist(tab(obj1))
tre <- njs(D)
jpeg(height=960, width=960)
par(xpd=TRUE, mar=c(0,0,0,0))
plot(tre, type="phylogram", edge.w=1)
dev.off()

# Plot observed vs expected heterozygosity for each locus
obj1sum <- summary(obj1)
jpeg(width=384, height=768)
plot(obj1sum$Hexp, obj1sum$Hobs, main="Observed vs expected heterozygosity by locus", xlab="Expected Heterozygosity", ylab="Observed Heterozygosity", xlim=c(0,0.5), ylim=c(0,1.0))
hmaxdf <- read.csv("HexpVsHmax.csv", header=TRUE) # a table of Hmax values for a range of Hexp values. I want to show what the maximum possible heterozygosity is for each expected heterozygosity.
lines(hmaxdf, type="l") #plots that limit line
lines(hmaxdf$Hmax, hmaxdf$Hmax, type="l", lty=2) #plot the line where Hobs = Hexp so it is easy to see which loci are above and below expected heterzygosity
dev.off()

#List the loci with more than 2 alleles
subset(obj1sum$loc.nall, obj1sum$loc.nall > 2)
# Output is 9 loci, each with 3 alleles (good, we expect few loci to have more than 2):
#L0046 L0072 L0424 L0607 L0749 L1042 L1479 L2089 L2188
#    3     3     3     3     3     3     3     3     3

#How many loci have zero observed heterozygosity?
length(subset(obj1sum$Hobs, obj1sum$Hobs == 0))
# Output is 756 loci
#[1] 756
#Get the locus names where Hobs = 0
H0names <- subset(names(obj1sum$Hobs), obj1sum$Hobs == 0)
#Write the names to a csv file so I can use it to filter the structure table in python
write.csv(H0names, "Hobs_0_locus_names.csv")

#How many loci have less than the expected H but not zero?
length(subset(obj1sum$Hobs, obj1sum$Hobs <= obj1sum$Hexp & obj1sum$Hobs > 0))
# Output
# [1] 97
# I can count the points that are above the Hobs=Hexp line but not on the Hmax line -
# there are only about 16 or 17 of them so about 115 points are on the Hmax line or Hobs=0 line.

#Doing this more precisely... find how many loci have 3 genotypes and how many have 2. Subtract the number with 0 heterozygosity from the number with 2 genotypes and this tells me exactly how many have Aa, AA only.
obj1df <- genind2df(obj1, sep="|") #make a genotypes dataframe
genotype_counts <- rapply(obj1df, function(x) length(unique(na.omit(x))))
length(subset(genotype_counts, genotype_counts == 3))
# Output
# [1] 119
length(subset(genotype_counts, genotype_counts == 2))
# Output
# [1] 2144
# (One locus, L0072 has 4 genotypes cause there are 3 different alleles)
# So... 2144 - 756 = 1388 that have the maximum possible heterozygosity (only Aa and AA)
#Get the list of loci with exactly 3 genotypes
G3names <- (names(subset(genotype_counts, genotype_counts == 3)))
#Write it to a csv file so I can use it to filter the structure table in python
write.csv(G3names, "Three_genotype_locus_names.csv")


#What is the distribution of observed Heterozygosity values?
summary(obj1sum$Hobs)
# Output:
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.00000 0.00000 0.03125 0.04060 0.03571 0.94120

#Write the genotypes dataframe to a csv file so I can use it to calculate MAF in python
write.csv(obj1df, file="Genotypes_dataframe_from_R.csv")

#Read the minor allele frequencies (generated from the file above using python) into a dataframe
mafs <- read.csv("minor_allele_frequencies.txt", header=FALSE)
summary(mafs)
# Output
#        V1
#  Min.   :0.01163
#  1st Qu.:0.01667
#  Median :0.02778
#  Mean   :0.03854
#  3rd Qu.:0.03448
#  Max.   :0.50000

#Plot histogram of MAF values
jpeg()
hist(mafs$V1, breaks=15, main="MAF distribution")
dev.off()


plot(Hetdf$Hexp, Hetdf$Hobs, pch=16, col=rgb(0, 0, 0, 0.035), main="Observed vs expected heterozygosity by locus", xlab="Expected Heterozygosity", ylab="Observed Heterozygosity", xlim=c(0,0.5), ylim=c(0,1.0))



#Import the structure table filtered (in python) to only include the loci with Hobs=0
obj1H0 <- read.structure("GT70_STRUCTURE_for_R_Hobs_0.str")
# It will prompt for info:
#   genotypes = 43 (number of samples)
#   markers = 756 (number of loci)
#   column with labels for genotypes = 1 (the names of the samples)
#   column with population factor = 0 (if we were specifying populations a priori but we aren't)
#   other optional columns - just hit enter, there are none in our data
#   row with the marker names = 0 (We don't have a marker names row)
#   genotypes coded by a single row = n (We have 2 lines per sample in the dataset)
indNames(obj1H0) # should return a list of the 43 sample names
ploidy(obj1H0) # should return 2 since we gave it 2 alleles for each marker

# Plot observed vs expected heterozygosity for each locus - just to confirm that they are all on the x-axis as I would expect if we've only included loci with no heterozygosity
obj1H0sum <- summary(obj1H0)
jpeg(width=384, height=768)
plot(obj1H0sum$Hexp, obj1H0sum$Hobs, main="Observed vs expected heterozygosity by locus", xlab="Expected Heterozygosity", ylab="Observed Heterozygosity", xlim=c(0,0.5), ylim=c(0,1.0))
hmaxdf <- read.csv("HexpVsHmax.csv", header=TRUE) # a table of Hmax values for a range of Hexp values. I want to show what the maximum possible heterozygosity is for each expected heterozygosity.
lines(hmaxdf, type="l") #plots that limit line
lines(hmaxdf$Hmax, hmaxdf$Hmax, type="l", lty=2) #plot the line where Hobs = Hexp so it is easy to see which loci are above and below expected heterzygosity
dev.off()

#List the loci with more than 2 alleles
subset(obj1H0sum$loc.nall, obj1H0sum$loc.nall > 2)
# Output is 2 loci, each with 3 alleles (good, we expect few loci to have more than 2):
# L354 L710
#   3    3

#Neighbor joining euclidian distance tree
DH0 <- dist(truenames(obj1H0))
treH0 <- nj(DH0)
jpeg(height=960, width=960)
par(xpd=TRUE, mar=c(0,0,0,0))
plot(treH0, type="phylogram", edge.w=1)
dev.off()

#Repeat for loci where Hobs = Hmax
obj1Hmax <- read.structure("GT70_STRUCTURE_for_R_Hmax.str")
# It will prompt for info:
#   genotypes = 43 (number of samples)
#   markers = 1391 (number of loci)
#   column with labels for genotypes = 1 (the names of the samples)
#   column with population factor = 0 (if we were specifying populations a priori but we aren't)
#   other optional columns - just hit enter, there are none in our data
#   row with the marker names = 0 (We don't have a marker names row)
#   genotypes coded by a single row = n (We have 2 lines per sample in the dataset)
indNames(obj1Hmax) # should return a list of the 43 sample names
ploidy(obj1Hmax) # should return 2 since we gave it 2 alleles for each marker
# Plot observed vs expected heterozygosity for each locus - just to confirm that they are all on the Hmax line
obj1Hmaxsum <- summary(obj1Hmax)
jpeg(width=384, height=768)
plot(obj1Hmaxsum$Hexp, obj1Hmaxsum$Hobs, main="Observed vs expected heterozygosity by locus", xlab="Expected Heterozygosity", ylab="Observed Heterozygosity", xlim=c(0,0.5), ylim=c(0,1.0))
hmaxdf <- read.csv("HexpVsHmax.csv", header=TRUE) # a table of Hmax values for a range of Hexp values. I want to show what the maximum possible heterozygosity is for each expected heterozygosity.
lines(hmaxdf, type="l") #plots that limit line
lines(hmaxdf$Hmax, hmaxdf$Hmax, type="l", lty=2) #plot the line where Hobs = Hexp so it is easy to see which loci are above and below expected heterzygosity
dev.off()
#Neighbor joining euclidian distance tree
DHmax <- dist(truenames(obj1Hmax))
treHmax <- nj(DHmax)
jpeg(height=960, width=960)
par(xpd=TRUE, mar=c(0,0,0,0))
plot(treHmax, type="phylogram", edge.w=1)
dev.off()


#Repeat for loci where all 3 genotypes are present
obj1G3 <- read.structure("GT70_STRUCTURE_for_R_3_genotypes.str")
# It will prompt for info:
#   genotypes = 43 (number of samples)
#   markers = 119 (number of loci)
#   column with labels for genotypes = 1 (the names of the samples)
#   column with population factor = 0 (if we were specifying populations a priori but we aren't)
#   other optional columns - just hit enter, there are none in our data
#   row with the marker names = 0 (We don't have a marker names row)
#   genotypes coded by a single row = n (We have 2 lines per sample in the dataset)
indNames(obj1G3) # should return a list of the 43 sample names
ploidy(obj1G3) # should return 2 since we gave it 2 alleles for each marker
# Plot observed vs expected heterozygosity for each locus - between Hobs=0 and Hmax
obj1G3sum <- summary(obj1G3)
jpeg(width=384, height=768)
plot(obj1G3sum$Hexp, obj1G3sum$Hobs, main="Observed vs expected heterozygosity by locus", xlab="Expected Heterozygosity", ylab="Observed Heterozygosity", xlim=c(0,0.5), ylim=c(0,1.0))
hmaxdf <- read.csv("HexpVsHmax.csv", header=TRUE) # a table of Hmax values for a range of Hexp values. I want to show what the maximum possible heterozygosity is for each expected heterozygosity.
lines(hmaxdf, type="l") #plots that limit line
lines(hmaxdf$Hmax, hmaxdf$Hmax, type="l", lty=2) #plot the line where Hobs = Hexp so it is easy to see which loci are above and below expected heterzygosity
dev.off()
#Neighbor joining euclidian distance tree
DG3 <- dist(truenames(obj1G3))
treG3 <- nj(DG3)
jpeg(height=960, width=960)
par(xpd=TRUE, mar=c(0,0,0,0))
plot(treG3, type="phylogram", edge.w=1)
dev.off()


# Import the table of geographic distances. I made this table in the same order as the data in the dist() table of genetic distances: S01 vs S02, S03....., S02 vs S03, S04......, S03 vs S04, S05..... and so on. Samples from the same collection (extractions from the same petri dish) have distance=0. Samples from different collections at the same site have distance = 0.2km. It is an arbitrary distance - enough to make them different from things in the same collection but much smaller than the smallest distance between two sites.
geog_dist <- read.csv("Crepidomanes_geographic_distance.csv", header=TRUE)
#print a summary of geographic distances
summary(geog_dist)
#Output:
# Distance_km
# Min.   :   0.0
# 1st Qu.:  15.6
# Median :  47.0
# Mean   : 401.6
# 3rd Qu.: 813.6
# Max.   :1504.7

#plot geographic distance vs genetic distance
geog_d <- geog_dist[1:903,]
genet_d <- D[1:903]
reg <- lm(genet_d~geog_d)
plot(geog_d, genet_d, main="Genetic vs geographic distance for each pair of samples", xlab="Geographic distance (km)", ylab="Genetic distance (Euclidian)")
abline(reg)

#Get r and r-squared values for the correlation
cor(geog_d, genet_d)
# Output (r):
# [1] -0.03703721
cor(geog_d, genet_d)^2
# Output (r-squared)
# [1] 0.001371755
#test the correlation
cor.test(geog_d, genet_d)
# Output
#    Pearson's product-moment correlation
#
# data:  geog_d and genet_d
# t = -1.1125, df = 901, p-value = 0.2662
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.10203001  0.02827044
# sample estimates:
#         cor
# -0.03703721

#Make a dist() object from my geographic distances. Easiest way is to make a copy of the genetic dist object and replace the values with the geographic distances
Dgen <- D
Dgeo <- D
Dgeo[1:903] <- geog_d

#Do a mantel test of isolation by distance
ibd <- mantel.randtest(Dgen, Dgeo)
ibd
#Output
# Monte-Carlo test
# Call: mantel.randtest(m1 = Dgen, m2 = Dgeo)
#
# Observation: -0.03703721
#
# Based on 999 replicates
# Simulated p-value: 0.596
# Alternative hypothesis: greater
#
#      Std.Obs  Expectation     Variance
# -0.322138451  0.001371325  0.014215766

#Plot the results of the mantel test
jpeg()
plot(ibd)
dev.off()

#Plot genetic distance vs geographic distance (with colors to indicate density of plotted points)
jpeg()
dens <- kde2d(Dgeo,Dgen, n=300, lims=c(-55, 1555,3.5,48))
myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
plot(Dgeo, Dgen, pch=20,cex=.5, xlab="Geographic distance (km)", ylab="Genetic distance")
image(dens, col=transp(myPal(300),.7), add=TRUE)
abline(lm(Dgen~Dgeo))
title("Isolation by distance plot")
dev.off()

#Repeat geographic vs genetic distance comparisons for just the 756 Hobs=0 loci
jpeg()
DgenH0 <- dist(truenames(obj1H0))
Dgeo <- DgenH0
Dgeo[1:903] <- geog_dist[1:903,]
myPal <- colorRampPalette(c("white","blue","gold", "red"))
dens <- kde2d(Dgeo,DgenH0, n=300, h=c(100,5), lims=c(-51, 1555,-1.2,36))
plot(Dgeo, DgenH0, pch=20,cex=.5, xlab="Geographic distance (km)", ylab="Genetic distance")
image(dens, col=transp(myPal(300),.7), add=TRUE)
abline(lm(DgenH0~Dgeo))
title("Isolation by distance plot")
dev.off()
ibd <- mantel.randtest(DgenH0, Dgeo)
ibd
# Output:
# Monte-Carlo test
# Call: mantel.randtest(m1 = DgenH0, m2 = Dgeo)
#
# Observation: -0.03944799
#
# Based on 999 replicates
# Simulated p-value: 0.601
#Alternative hypothesis: greater
#
#      Std.Obs  Expectation     Variance
# -0.354374480  0.003352258  0.014587055
jpeg()
plot(ibd)
dev.off()

#Get IBD plots for each sample individually (i.e., for each sample plot it's genetic vs geographic distance to each other sample). This may help me spot samples that have a very different pattern from everything else.
for (i in 1:43) { plot(Dgeo_df[,i], DgenH0_df[,i], xlim=c(0,1550), ylim=c(0,35), xlab="Geographic distance", ylab="Genetic Distance", main=names(DgenH0_df)[i]) }


#import the structure table filtered for only Hobs=0 and samples with low (<11) average genetic distance
obj1H0_low <- read.structure("GT70_STRUCTURE_for_R_Hobs_0-lowDgen.str")
# It will prompt for info:
#   genotypes = 26 (number of samples)
#   markers = 756 (number of loci)
#   column with labels for genotypes = 1 (the names of the samples)
#   column with population factor = 0 (if we were specifying populations a priori but we aren't)
#   other optional columns - just hit enter, there are none in our data
#   row with the marker names = 0 (We don't have a marker names row)
#   genotypes coded by a single row = n (We have 2 lines per sample in the dataset)
indNames(obj1H0_low) # should return a list of the 26 sample names
ploidy(obj1H0_low) # should return 2 since we gave it 2 alleles for each marker

DH0_low <- dist(truenames(obj1H0_low))
treH0_low <- nj(DH0_low)
par(xpd=TRUE, mar=c(0,0,0,0))
plot(treH0_low, type="phylogram", edge.w=1)

DgenH0_low <- dist(truenames(obj1H0_low))
Dgeo_low <- DgenH0_low
geog_dist_low <- read.csv("Crepidomanes_geographic_distance-LowDgen.csv", header=TRUE)
Dgeo_low[1:325] <- geog_dist_low[1:325,]
myPal <- colorRampPalette(c("white","blue","gold", "red"))
dens <- kde2d(Dgeo_low,DgenH0_low, n=300, h=c(100,5), lims=c(-51, 1555,-1.2,36))
plot(Dgeo_low, DgenH0_low, pch=20,cex=.5, xlab="Geographic distance (km)", ylab="Genetic distance")
image(dens, col=transp(myPal(300),.7), add=TRUE)
abline(lm(DgenH0_low~Dgeo_low))
title("Isolation by distance plot")
ibd <- mantel.randtest(DgenH0_low, Dgeo_low)
ibd
#Output
# Monte-Carlo test
# Call: mantel.randtest(m1 = DgenH0_low, m2 = Dgeo_low)
#
# Observation: -0.04916927
#
# Based on 999 replicates
# Simulated p-value: 0.639
# Alternative hypothesis: greater
#
#      Std.Obs  Expectation     Variance
# -0.407448337 -0.002636675  0.013042760
plot(ibd)
