###set working directory

#install and load packages
library(ggplot2)
library(viridis)

#import the data
#This can be changedfor the dataset
df<- read.table('dros_comparison_high.txt',sep='\t',header=T)

#This creates the plot, the first line loads in the data
ggplot(df,aes(fill=spec, y=val, x=file))+
  #states that it is a bar plot
  geom_bar(position = 'dodge',stat='identity')+
  #adds labels
  labs(ylab = "CpG dyads per 1000 BP", title='A', size= 200)+
  #sets colours
  scale_fill_viridis_d() +
  #more labs
  labs(x = "Mouse Gene", y = "CpG codon dyad frequency")+
  #general themeing
  theme(legend.position = "none",  panel.background = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(color = "black", size = 16), text=element_text(size=20), plot.title=element_text(color = "black", size = 40))
  

#A t-test across the data 
overall_t_test <- t.test(df$val[df$spec == "dros"], df$val[df$spec == "mal"],paired=T)

overall_t_test

#wilcox test across the data
overall_w_test <- wilcox.test(df$val[df$spec == "dros"], df$val[df$spec == "mal"], paired=T)
overall_w_test