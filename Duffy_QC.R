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
#   genotypes = 244  (number of samples) This number can be found in the ipyrad _stats file, I had 266 but I threw out two samples with no data. 
#   markers = 1886 (number of loci) Also find in ipyrad _stats file.
#   column with labels for genotypes = 1 (the names of the samples)
#   column with population factor = 0 (if we were specifying populations a priori but we aren't)
#   other optional columns - just hit enter, there are none in our data
#   row with the marker names = 0 (We don't have a marker names row)
#   genotypes coded by a single row = n (We have 2 lines per sample in the dataset.)
indNames(obj1) # should return a list of the 244 sample names, I only get 229 because low coverage individuals get filtered out.
ploidy(obj1) # should return 2 since we gave it 2 alleles for each marker.


#Neighbor joining euclidian distance tree
D <- dist(tab(obj1))               #super hard to read, create a distance matrix! 

#What does D look like?
D                                  #0 means they are identical...and should be expected across the diagonal as each sample is identical to itself. 0s other than the diagonal better be replicates... Nas just mean that there are no loci in common between the two samples (we have no info ont their relatedness.)
class(D)                           #dist
attributes(D)
str(D)
D[100]                             #can index by row number
D["p_019s_11"]                     #can index by sample name

#Convert to matrix for indexing.
M <- as(D, "matrix")               #Makes a normal matrix of the distance values.
M                                  
str(M)

#Indexing Distance Matrix (D)
#All of the distance values for each individual are in each row, and each column, due to nature of the matrix.
row1 <- M["p_001s_01",]
row1                                #row 1 has the distance values between p_001s_01 and every other sample in this library.
str(row1)
row1["p_001s_03"]                   #index by sample name
row1[3]                             #by index number

#Function for identification of the 5 closest relatives to a sample.
  #What I need to do:
    # in order to find the 3 lowest, find the lowest, set it to an infinitely high value, then find the lowest again, repeat
  #pseudocode:
    #for "p001-s_01" in 1 to 229 in row
    #index number = decreasing_index[i]
    #row1[index number] 
    #if NA, we stop
    #if first time through, save, its the lowest value!
    #as long as it is 0, add index number to the list
    #then we have a list that is second lowest
    #output will be a list of index numbers, then search names_list to get the sample name for that index

#Start small, write the code for just row 1:
find_relatives <-function(row){
  relatives <- list()                                                        #initialize an empty list
  decreasing_index <- order(row, decreasing = FALSE)                         #sort, smallest distance values first
  for (i in 1:5){
    relatives[i] <- decreasing_index[i]                                      #first item in the list is the closest relative
  }
  return(relatives)
}
relatives <- find_relatives(row1)

#Now I want the output of find_relatives to be converted into sample names, fount in names_list.
index_to_samples <- function(relatives){
  samples <- list()
  for(i in relatives){
    print(names_list[i])
  }
}
index_to_samples(relatives)



#Possibly helpful:
which.max(row1)
order(row1)
decreasing_index <- order(row1, decreasing = FALSE)
decreasing_index[229]
which(row1 == min(row1))
my_list = c(1,2,3)
my_list[2]

#Making a Tree
tre <- njs(D)
jpeg(height=960, width=960)
par(xpd=TRUE, mar=c(0,0,0,0))
plot(tre, type="phylogram", edge.w=1)
dev.off()