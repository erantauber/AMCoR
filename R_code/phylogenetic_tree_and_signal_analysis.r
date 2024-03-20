###set your working directory
#install and load packages
library(ape)
library(phytools)
library(geiger)
library(phylobase)
library(picante)
library(adephylo)
library(ade4)

#load in the data
tree <- read.nexus('261Mammals.consensus.tree')
data <- read.csv('species_list_with_data.txt',sep='\t')




#Creates a cols scale that is used in the graph
cols2 <- viridis::inferno(length(unique(rounded)+20))
#an empty colour vector which will contain each tips colour code
cols=c()
#for loop that checks the value of each tip label and colours it accordingly.
for (i in data$rounded){
  my_num <- (i/10000)-24
  temp_col<- cols2[my_num]
  cols<-c(cols,temp_col)
}
#re-rooting of the tree to the oldest species
tree_root <- 'Solenodon_paradoxus'
tree <- root(tree, tree_root)
#converting the tree to an ultrametric one (optional)
ultrametric_tree <- chronos(tree)
#ploting the tree. 
plot(ultrametric_tree, tip.color = cols,cex=.6, type = 'fan')

###below is the starts into the model type of the data, as well as our investigation into phylogenetic signal. 

#This creates an vector that is used by the fitcontinuous below
nums<-setNames(data$num,tree$tip.label)

#This is the three models under our analysis
fitBM<-fitContinuous(tree,nums)
fitOU<-fitContinuous(tree,nums,model="OU")
fitEB<-fitContinuous(tree,nums,model="EB")

#This pulls the AIC values from our models
aic.vals<-setNames(c(fitBM$opt$aicc,fitOU$opt$aicc,fitEB$opt$aicc),
                   c("BM","OU","EB"))
#This prints off the results of the analysis
fitBM
fitOU
fitEB
#This prints the Aic values for our different models for comparison
aic.vals
#This weighs up the different models to determine which is the most accurate of those tested
aic.w(aic.vals)

#The tests below are for phylogenetic signal

phylotraits <- phylo4d(tree, data)
moran.test <- abouheif.moran(phylotraits,method="Abouheif")
moran.test

#Abouheif model
abouheif.test <- abouheif.moran(phylotraits,method="oriAbouheif")
abouheif.test


#lamda model
phylosig(tree, nums, method="lambda", test=TRUE, nsim=999)