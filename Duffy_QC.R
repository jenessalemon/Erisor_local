setwd('/Users/jimblotter/Desktop/Grad_School/Data_Analysis/QC/') #remember to put the .str file there!
#install.packages("genetics")
library("ape")
library("genetics")        #there is no package called ‘genetics’ = install.packages("genetics")
library("pegas")
library("seqinr")
library("ggplot2")
library("adegenet")

#Read in data
obj1 <- read.structure("eri_sor_9sout.str", n.ind = 244, n.loc = 1886, col.lab = 1, col.pop = 0, col.others = NULL, row.marknames = 0) #place cursor in console
# It will prompt for info:
#   genotypes = 246 244 because I threw out two (number of samples) GOT THIS FROM BARCODES FILE
#   markers = 1886 (number of loci)
#   column with labels for genotypes = 1 (the names of the samples)
#   column with population factor = 0 (if we were specifying populations a priori but we aren't)
#   other optional columns - just hit enter, there are none in our data
#   row with the marker names = 0 (We don't have a marker names row)
#   genotypes coded by a single row = n (We have 2 lines per sample in the dataset? I THINK I JUST HAVE ONE)
indNames(obj1) # should return a list of the 245 sample names, I only get 229
ploidy(obj1) # should return 2 since we gave it 2 alleles for each marker DOES

#Neighbor joining euclidian distance tree
D <- dist(tab(obj1))

#What does D look like?
D
dim(D)                             #NULL, Apparenty because this is a "dist" object.
class(D)                           #dist
attributes(D)
str(D)
D[1]                               #didnt get a 0 so that is good
which(D==62.34581)                 #didn't work. integer(0)

attributes(dist1)
M <- as(D, "matrix")               #this maybe makes a normal matrix of the distance values?
M                                  #there's a lot of 0s, but there are actually values appearing?

#Indexing Distance Matrix (D)
D

#Making a Tree
tre <- njs(D)
jpeg(height=960, width=960)
par(xpd=TRUE, mar=c(0,0,0,0))
plot(tre, type="phylogram", edge.w=1)
dev.off()