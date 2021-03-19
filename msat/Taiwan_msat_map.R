---
  title: "Taiwan Kudzu pop gen"
author: "Sandra Hoffberg"
date: "Oct 31,2019"
output: html_document
---
  
  
  ```{r load package}
#first load libraries
library(poppr)
devtools::install_github("dwinter/mmod", ref="bs_fix")
library(mmod)

#import microsat dataset that has GPS coordinates
Taiwanmsats<-read.genalex("C:/Users/Sandra/OneDrive/kudzu/Taiwan/Taiwan_msat_GPS_least_missing2.csv", 
                          geo=T, region=F, ploidy=2, recode=F, genclone=T)

#look at data
Taiwanmsats

```{r select clonal threshold}
# Repeat lengths are necessary for Bruvo's distance
msatsize<-c(4,4,4,4,4,4,3,4,5,5,4,3,3)

# Now we can collect information of the thresholds. We can set threshold = 1 because we know that this will capture
# the maximum possible distance:
invisible(thresholds <- mlg.filter(Taiwanmsats, distance = bruvo.dist, stats = "THRESHOLDS", replen = msatsize, 
                                   threshold = 1, algorithm = "average_neighbor", missing="ignore"))

# We can use these thresholds to find an appropriate cutoff
pcut <- cutoff_predictor(thresholds)
pcut
{hist(thresholds, nclass=200)
  abline(v=pcut, col="red")
  abline(v=.01, col="blue")}

mlg.filter(Taiwanmsats, distance = bruvo.dist, replen = msatsize, algorithm="average_neighbor", missing="ignore") <- 0.01
Taiwanmsats

#At this point clones are defined, but the genotypes are not not changed.  so the dataset still has whatever somatic mutations and errors in it.  It still has the same number of samples, so it is NOT a genet-only dataset.   

Taiwanmsat2<-clonecorrect(Taiwanmsats, strat=~Pop, combine = F, keep = 1)
Taiwanmsat2

```{r DAPC + map}
library(adegenet)
library(ape)
library(ggplot2)
library(devtools)
library(dplyr)
library(stringr)
library(maps)
library(mapdata)
library(mapplots)
library(scatterpie)



#CHoose number of clusters 
#compare between clone corrected and full dataset

grp<-find.clusters(Taiwanmsats) #keep all PCs
150
5

#Keep all PCS and then keep # clusters (elbow)

grp
table(pop(Taiwanmsats), grp$grp)
table.value(table(pop(Taiwanmsats), grp$grp), col.lab=paste("inf", 1:50),
            row.lab=paste("ori", 1:45))


grp2<-find.clusters(Taiwanmsat2) #keep all PCs
150
5

grp2
table(pop(Taiwanmsat2), grp2$grp)
table.value(table(pop(Taiwanmsat2), grp2$grp), col.lab=paste("inf", 1:6),
            row.lab=paste("ori", 1:45))


#Run DAPC

dapc.Taiwan <- dapc(Taiwanmsats, grp$grp)
40
15
#Do not keep all PCs. Can keep all in barplot if small #
scatter(dapc.Taiwan, scree.da=FALSE, bg="white", pch=20, cell=0, cstar=0, col=rainbow(5), solid=.4,
        cex=3,clab=0, leg=TRUE, txt.leg=paste("Cluster",1:5))

something <- optim.a.score(dapc.Taiwan)
# 17 and up is good


dapc.Taiwan2 <- dapc(Taiwanmsat2, grp2$grp)
40
4
#Do not keep all PCs. Can keep all in barplot if small #
scatter(dapc.Taiwan2, scree.da=FALSE, bg="white", pch=20, cell=0, cstar=0, col=rainbow(5), solid=.4,
        cex=3,clab=0, leg=TRUE, txt.leg=paste("Cluster",1:5))


summary(dapc.Taiwan2)


something <- optim.a.score(dapc.Taiwan2)
# ~7 is good


pramx <- xvalDapc(tab(Taiwanmsats, NA.method = "mean"), pop(Taiwanmsats))
system.time(pramx <- xvalDapc(tab(Taiwanmsats, NA.method = "mean"), #pop(Taiwanmsats),
                              grp=grp$grp,
                              n.da=5,
                              #n.pca = 6:14,
                              n.rep = 100,
                              parallel = "multicore", ncpus = 4L))
pramx

set.seed(999)

pramx2 <- xvalDapc(tab(Taiwanmsat2, NA.method = "mean"), pop(Taiwanmsat2))
system.time(pramx2 <- xvalDapc(tab(Taiwanmsat2, NA.method = "mean"), #pop(Taiwanmsat2),
                               grp=grp2$grp,
                               n.da=5,
                              #n.pca = 6:14,
                              n.rep = 100,
                              parallel = "multicore", ncpus = 4L))

pramx2




scatter(pramx$DAPC, cstar=0)
scatter(pramx2$DAPC, cstar=0)

pramx2$DAPC$posterior

#col=c("blue", "black", "forestgreen", "red", "gray", "magenta", "gold", "skyblue", "purple")

compoplot(pramx2$DAPC)
```


######## R map with clones ######
world <- map_data("world")
Taiwan<-subset(world, region %in% c("China", "Taiwan"))


Taiwan<-ggplot() + 
  geom_polygon(data = Taiwan, aes(x=long, y = lat, group = group), fill = NA, color = "black") + 
  coord_map(xlim = c(100,123),ylim = c(18,30)) # coord_fixed(1.3)


GPS_Taiwan<-read.table("C:/Users/Sandra/OneDrive/kudzu/Taiwan/GPS_by_pop_msat.txt", header = T)
GPS_Taiwan
length(GPS_Taiwan[,1])

labs <- data.frame(
  long=GPS_Taiwan[,3],
  lat=GPS_Taiwan[,2],
  names=GPS_Taiwan[,1],
  stringsAsFactors = FALSE
); labs

data <-as.data.frame(unclass(table(pop(Taiwanmsats), grp$grp))); data
colnames(data)<- c("Group 1","Group 2","Group 3","Group 4","Group 5")#,"Group 6","Group 7",
                   #"Group 8","Group 9","Group 10")#,"Group 11","Group 12","Group 13","Group 14",
                   #"Group 15")


colors3<-rainbow(5)#, s = 1, v = 1, start = 0, end = max(1, 20 - 1)/20, alpha = 1)#, rev = FALSE)


Taiwan_data<-cbind(GPS_Taiwan[,2:3], data); Taiwan_data

#put sampling points on the map
Taiwan + 
  geom_point(data = labs, aes(x = long, y = lat), color = "black", size =2) 



Taiwan +
  geom_scatterpie(data=Taiwan_data, aes(x=Long, y=Lat, r=.2),sorted_by_radius = F, cols= colnames(data)) +
  scale_fill_manual(values=c("red","hotpink","blue","orange","maroon4")) #rainbow(8))

ggsave("C:/Users/Sandra/OneDrive/kudzu/Taiwan/Taiwan_msat_map_withclones.svg", device="svg", dpi=600)


######## R map without clones ######

world <- map_data("world")
Taiwan<-subset(world, region %in% c("China", "Taiwan"))


Taiwan<-ggplot() + 
  geom_polygon(data = Taiwan, aes(x=long, y = lat, group = group), fill = NA, color = "black") + 
  coord_map(xlim = c(100,123),ylim = c(18,30)) # coord_fixed(1.3)


GPS_Taiwan<-read.table("C:/Users/Sandra/OneDrive/kudzu/Taiwan/GPS_by_pop_msat.txt", header = T)
GPS_Taiwan
length(GPS_Taiwan[,1])

labs <- data.frame(
  long=GPS_Taiwan[,3],
  lat=GPS_Taiwan[,2],
  names=GPS_Taiwan[,1],
  stringsAsFactors = FALSE
); labs

data <-as.data.frame(unclass(table(pop(Taiwanmsat2), grp2$grp))); data
colnames(data)<- c("Group 1","Group 2","Group 3","Group 4","Group 5")#,"Group 6","Group 7",
#"Group 8","Group 9","Group 10")#,"Group 11","Group 12","Group 13","Group 14",
#"Group 15")


colors3<-rainbow(5)#, s = 1, v = 1, start = 0, end = max(1, 20 - 1)/20, alpha = 1)#, rev = FALSE)


Taiwan_data<-cbind(GPS_Taiwan[,2:3], data); Taiwan_data

#put sampling points on the map
Taiwan + 
  geom_point(data = labs, aes(x = long, y = lat), color = "black", size =2) 



Taiwan +
  geom_scatterpie(data=Taiwan_data, aes(x=Long, y=Lat, r=.4),sorted_by_radius = F, cols= colnames(data)) +
  scale_fill_manual(values=c("orange","blue","maroon4","hotpink","red")) #rainbow(8))

ggsave("C:/Users/Sandra/OneDrive/kudzu/Taiwan/Taiwan_msat_map_without_clones.svg", device="svg", dpi=600)
