library(ape)
library(apTreeshape)
library(adephylo)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(broom)
library(zoom)

#
#set working directory
setwd("~/Desktop/exam")

#read tree file
Monkeytree<- read.nexus("alignmaf_3.nexus.con.tre")

#import time
time <- read.table("monkeypox_db.csv", header = TRUE, sep = ",", nrows = 32, quote="",
                   colClasses = c("NULL", "NULL", "NULL", "NULL", "numeric"))
names(time) <- c('Collection_Date')
#find distance from root (using adephylo package)
distancefromroot= distRoot(Monkeytree, tips = "all", method = c("patristic"))

#find distance from root (using ape)
distance_all=cophenetic.phylo(Monkeytree)
distancefromoutgroup=distance_all[,48]
apedistancehumans=distancefromoutgroup[1:37]


#create dataframe with only human hosts
df <- data.frame(Time =  time,
                 Distance =  apedistancehumans
                 )
df2 <- data.frame(Distance = covidtree$edge.length)
corplot=ggplot(data = df, aes(x = Time, y = Distance)) +
geom_point()
print(corplot + labs(y = "Distance", x = "Time (days)"))

#linear regression analysis
linearr <- lm(Distance ~ Time, data = df)
#extract coefficients
coef1=linearr$coefficients[1]
coef2=linearr$coefficients[2]
#take a printscreen
summary(linearr)
s=summary(linearr)
#draw the line
y=linearr$coefficients[1]+linearr$coefficients[2]*df$Time
df$model=y

corplotlr=ggplot(data = df, aes(x = Time, y = Distance)) +
  geom_point() + geom_line(aes(y = model), color = "red", linetype = "dotted") +
                             geom_text(    label=rownames(df),      check_overlap = T , nudge_x = 10,nudge_y=0.00005)
  
print(corplotlr + labs(y = "Distance", x = "Time (days)"))

#plot tree and find the distance of last human infected ancestor
plot(Monkeytree, show.tip.label=TRUE)
add.scale.bar(length=0.005)
edgelabels(cex=0.8, frame = 'none')
zm()
#distance of last human COVID19 ancestor is equal to total distance - 
lastnode=sum(df2$Distance[1:6])
#0.45863=4.587e-01+2.502e-06*df$Time
#-4.587e-01=2.502e-06*df$Time
#-4.587e-01/2.502e-06=TimeBefore
Initialmutation=(lastnode-coef1)/coef2
Initialmutation
