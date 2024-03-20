###set working directory

#Edit this to change which data files you are processing 
file<-'NNCGNN'

####Supplementary graph one 
#The file extension is combined with the file data above to open the correct file
extension<-'_conserved_locations.csv'
#This table provides the data for the red line that will overlay the graph indicating the percentage of genes that are at least that length
survival <- read.table('gene_length_survival.txt',sep='\t', header=T)

#double check that the data is imported, or you will just be using the old data. 
df1<- read.csv(paste(file,extension,sep=''))
#As there are very few data points for genes over 8000 bp long, these are removed from the graph to improve visualisation. This cut off can be altered if you would like.
my_ind<- which(df1$location>=8000)[1]
df1<-df1[c(0:my_ind),]
#The code not only makes the graph, but the tiff file at the end. To display the graph in R, skip the 'tiff()' line
file_name<-paste(file,'_frequency_of_CpG.tiff')
#If you want to create a file automatically, use this. Otherwise skip so it is only displayed in R
tiff(file_name,width=2000,height = 1600,units = "px", res = 300)
#This creates the barplot 
barplot(df1$no.of.conserved.dyads~df1$location,
        ylab='Frequency',
        xlab='Dyad Position')
#As we are going to overlay a new graph we need to set this up with par()
par(new = TRUE)
#Again, the data needs to be trimmed to match the previous dataset
survival<- survival[c(0:my_ind),]
#This plots a line graph on top of the previous one. Crucially the axis is not added here, as we want it to be on the right.
plot(survival$loc,survival$survivors, type = "l",  pch = 16, 
     xlab = "", ylab = "", axes = FALSE, col='red')
#Axis is added
axis(4, at = seq(0, max(survival$survivors), by = 10), pos= 8000,
     labels = seq(0, max(survival$survivors ), by = 10))
#Axis lab is added. 
mtext("Percentage of Genes", side = 4, line = 0.9)
#This ends the creation of the plot, and creates the image file stated in the tiff() line
dev.off()


#Figure 2

#Data is loaded in
extension<- '_percentage_locations.csv'
df2<- read.csv(paste(file,extension,sep = ''))

#Due to the way that python handle floating points, we sometimes get percentages above 100. These are trimmed.
df2<-df2[c(0:100),]
#creates the general percentage data for all genes
file_name<-paste(file,'_percetage_location_graph.tiff')
#Same process as the above for generating the graph as well as the tiff file
tiff(file_name,width=2000,height = 1600,units = "px", res = 300)
barplot(df2$Freq~df2$Percent,
        xlab = 'Percentage location',
        ylab = 'Frequency',  cex.lab = 1.5)
dev.off()


#supplementary figure 2

#This follows much the same process as above
df3 <- read.csv(paste(file,'_top_1000_percentage_locations.csv',sep=''))
df3<-df3[c(0:100),]
file_name<-paste(file,'_top_1000_percetage_location_graph.tiff')
tiff(file_name,width=2000,height = 1600,units = "px", res = 300)
barplot(df3$Freq~df3$Percent,
        xlab = 'Percentage location',
        ylab = 'Frequency',  cex.lab = 1.5)
dev.off()



df3 <- read.csv(paste(file,'_bottom_1000_percentage_locations.csv',sep=''))
df3<-df3[c(0:100),]
file_name<-paste(file,'_bottom_1000_percetage_location_graph.tiff')
tiff(file_name,width=2000,height = 1600,units = "px", res = 300)
barplot(df3$Freq~df3$Percent,
        xlab = 'Percentage location',
        ylab = 'Frequency',  cex.lab = 1.5)
dev.off()