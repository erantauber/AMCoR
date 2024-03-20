###set working directory

#load in data
df<-read.table('mice_stage_data.txt', sep = '\t')

#this generated the basic plot
plot(df$V2,
     xlab = 'Developmental Stage',
     ylab= 'CpG Dyad Value', type = 'p', pch=16,
     ylim=c(0.16,0.21),
     xlim=c(.5,9))

#This adds the names of the stages, offset by a small margin for clarity
for (nam in row(df)){
  if (nam==4 || nam==2){
    text(nam,df[nam,2]-0.003,df[nam,1], cex=0.9)
  }
  else{
  text(nam,df[nam,2]+0.003,df[nam,1], cex=0.9)
  }
  if (nam==8){
    break
  }
}
#This creates the abline that is printed over the top. 
abline(lm(df$V2~c(1:8)), col='red')

