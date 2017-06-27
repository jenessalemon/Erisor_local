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

names_list <- as.list(indNames(obj1))
names_list[4]


#Neighbor joining euclidian distance tree
D <- dist(tab(obj1))

#What does D look like?
D                                  #this is strange looking. Is it ok that there are NA and 0s? 0 means they are identical... and there are so many blank spaces...?
dim(D)                             #NULL, Apparenty because this is a "dist" object.
class(D)                           #dist
attributes(D)
str(D)
D[100]                               #didnt get a 0 so that is good
D["p_019s_11"]
D.labels

attributes(dist1)
M <- as(D, "matrix")               #this maybe makes a normal matrix of the distance values?
M                                  #there's a lot of 0s, but there are actually values appearing?

#Indexing Distance Matrix (D)
  #pseudocode:
  #for a certain individual (replicate):
    # for each row, find the lowest 3 nubmers, return the index for that item
    # for the 2nd column, we need everything before the second item, and everything after.
    # in order to find the 3 lowest, find the lowest, remove it, then find the lowest again, repeat

row1 <- M["p_001s_01",]
row1
row1["p_001s_03"]
str(row1)
which.max(row1)
row1[226]
order(row1)
decreasing_index <- order(row1, decreasing = FALSE)
row1[1]
which(row1 == min(row1))

#for "p001-s_01" in 1 to 229 in row
#index number = decreasing_index[i]
#row1[index number] 
#if NA, we stop
#if first time through, save, its the lowest value!
#as long as it is 0, add index number to the list
#then we have a list that is second lowest

#output will be a list of index numbers, then search names_list to get the sample name for that index

#Making a Tree
tre <- njs(D)
jpeg(height=960, width=960)
par(xpd=TRUE, mar=c(0,0,0,0))
plot(tre, type="phylogram", edge.w=1)
dev.off()