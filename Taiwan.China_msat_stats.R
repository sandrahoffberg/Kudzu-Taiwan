#first load libraries
library(poppr)
devtools::install_github("dwinter/mmod", ref="bs_fix")
library(mmod)


#import microsat dataset that has GPS coordinates
Taiwanmsats<-read.genalex("C:/Users/Sandra/OneDrive/kudzu/Taiwan/Taiwan_msat_GPS_least_missing2_noKYU.csv", geo=T, region=F, ploidy=2, recode=F, genclone=T)

#look at data
Taiwanmsats
summary(Taiwanmsats)


# Do I have enough loci to distinguish individuals? 
set.seed(2102021)
genotype_curve(Taiwanmsats, sample=1000, maxloci=13, thresh = 0.9)
# getthelastplot
p <-last_plot()+theme_bw()
svg("C:/Users/Sandra/OneDrive/kudzu/Taiwan/msat_info.svg")
p
dev.off()
# Yes


#check to see if we need to do a clonal correction
#Psex is the probability of encountering a genotype more than once by chance
#If Psex < 0.05, duplicated genotypes should be treated as different individuals. In contrast, if Psex is > 0.05, then you would treat the two genotypes as clones.
psex1<-psex(Taiwanmsats, by_pop=T,method="multiple")

#view the results by color with probability of matching samples being clones in red
plot(psex1, col = ifelse(psex1 > 0.05, "red", "black"), pch=20)

#find the maximum of Psex. If it is higher than 0.05, do a clonal correction.
max(psex1)
# A clonal correction is needed. 


# Condense Kudzu dataset into clonal lineages in poppr using mlg.filter

#Thresholds for collapsing multilocus genotypes should not be arbitrary.
#The Poppr method of choosing a threshold is to find a gap in the distance distribution that represents clonal groups.
#You can look at this by analyzing the distribution of all possible thresholds with the function "cutoff_predictor".

#Start with a distance matrix.  Bruvo distance is made for msats, but it calculates the average over all loci in the population.
#This doesn't give you the threshhold necessary for mlg.filter. Therefore, we have to use cutoff_predictor. 


#Then you can use clonecorrect to actually remove repeated genotypes, and get a genet-only dataset, BUT I DID NOT DO THIS.


# Repeat lengths are necessary for Bruvo's distance
msatsize<-c(4,4,4,4,4,4,3,4,5,5,4,3,3)

# Now we can collect information of the thresholds. We can set threshold = 1 because we know that this will capture the maximum possible distance:
invisible(thresholds <- mlg.filter(Taiwanmsats, distance = bruvo.dist, stats = "THRESHOLDS", replen = msatsize, threshold = 1, algorithm = "average_neighbor", missing="ignore"))

# We can use these thresholds to find an appropriate cutoff
pcut <- cutoff_predictor(thresholds)
pcut

svg("C:/Users/Sandra/OneDrive/kudzu/Taiwan/msat_clonal_threshold.svg")
{hist(thresholds, nclass=200)
  abline(v=pcut, col="red")}
#  abline(v=.01, col="blue")}
dev.off()

# assign clones 
mlg.filter(Taiwanmsats, distance = bruvo.dist, replen = msatsize, algorithm="average_neighbor", missing="ignore") <- pcut
Taiwanmsats

#At this point clones are defined, but the genotypes are not not changed.  
#so the dataset still has whatever somatic mutations and errors in it.  
#It still has the same number of samples, so it is NOT a genet-only dataset.   




## Next, we get summary stats for each population. The output is: 

### A data frame with populations in rows and the following columns:
Pop - A vector indicating the population factor
N - An integer vector indicating the number of individuals/isolates in the specified population.
MLG - An integer vector indicating the number of multilocus genotypes found in the specified population, (see: mlg)
eMLG - The expected number of MLG at the lowest common sample size (set by the parameter minsamp).
SE - The standard error for the rarefaction analysis
H - Shannon-Weiner Diversity index
G - Stoddard and Taylor's Index
lambda - Simpson's index
E.5 - Evenness
Hexp - Nei's gene diversity (expected heterozygosity)
Ia - A numeric vector giving the value of the Index of Association for each population factor, (see ia).
p.Ia - A numeric vector indicating the p-value for Ia from the number of reshufflings indicated in sample. Lowest value is 1/n where n is the number of observed values.
rbarD - A numeric vector giving the value of the Standardized Index of Association for each population factor, (see ia).
p.rD - A numeric vector indicating the p-value for rbarD from the number of reshuffles indicated in sample. Lowest value is 1/n where n is the number of observed values.
File - A vector indicating the name of the original data file.'

```{r summary stats}
########### Genetic Diversity Stats #########
#MLG, eMLG +/- SE, Hexp +/- SE

summarystatT<-poppr(Taiwanmsats)
summarystatT
#write.csv(summarystat, "D:/kudzu2019/msat/Taiwan_msats_summary_stats.csv")
```

```{r Fst}
diff_stats(Taiwanmsats)
Fst<-pairwise_Gst_Nei(Taiwanmsats, linearized = FALSE) # Calculates pairwise Gst. If linearized = TRUE, it calculates 1/(1- Gst)
Fst2<-pairwise_Gst_Nei(Taiwanmsats, linearized = T) # Calculates pairwise Gst. If linearized = TRUE, it calculates 1/(1- Gst)  
Fst
Fst2

write.csv(as.matrix(Fst), "C:/Users/Sandra/OneDrive/kudzu/Taiwan/Tw.Ch_Gst_msat.csv")
write.csv(as.matrix(Fst2), "D:/kudzu2019/Taiwan/1-Gst_msat.csv")

bs <- chao_bootstrap(Taiwanmsats, nreps = 100)
summarise_bootstrap(bs, Gst_Nei)     # for Nei's Gst
```



#AMOVA
Taiwanmsats
Taiwanmsat2<-clonecorrect(Taiwanmsats, strat=~Pop, combine = F, keep = 1)
distmat2<-bruvo.dist(Taiwanmsat2, msatsize, add=TRUE, loss=TRUE, by_locus = F)

AMOVA2<- poppr.amova(Taiwanmsat2, hier = ~Pop, clonecorrect = T, within = F,
                     #dist = as.matrix(distmat2),  squared = F, freq = TRUE,
                     #correction = "quasieuclid", filter = FALSE,
                     threshold = .01, algorithm = "average_neighbor", threads = 2L,
                     missing = "loci", cutoff = 0.15, quiet = FALSE, method = "ade4", nperm = 0)
AMOVA2


distmat1<-bruvo.dist(Taiwanmsats, msatsize, add=TRUE, loss=TRUE, by_locus = F)
as.matrix(distmat1)
AMOVA1<- poppr.amova(Taiwanmsats, hier = ~Pop, clonecorrect = F, within = F,
                     #dist = as.matrix(distmat1),  squared = F, #freq = TRUE,
                     #correction = "quasieuclid", filter = FALSE,
                     threshold = 0, algorithm = "average_neighbor", threads = 2L,
                     missing = "loci", cutoff = 0.15, quiet = FALSE, method = "ade4", nperm = 0)

AMOVA1
```






#Now let's make some trees and start to figure out the invasion history!
  ```{r load libs}
library(phangorn)
````


```{r nj tree}
#Make a distance matrix using bruvo.dist
distmat<-bruvo.dist(Taiwanmsats, msatsize, add=TRUE, loss=TRUE, by_locus = F)
distmat
# View each population as a heatmap.
sapply(popNames(Taiwanmsats), function(x)
  heatmap(as.matrix(bruvo.dist(popsub(Taiwanmsats, x), replen = msatsize)), symm=TRUE))


njtree<-nj(distmat)
njtree
write.tree(njtree, "D:/kudzu2019/Taiwan/NJ_msat_Taiwan.tre")

njnet<-neighborNet(distmat)
```







```{r DAPC + map}
##########see R script for this, this was re-done


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

grp<-find.clusters(Taiwanmsat2) #keep all PCs
#40 then 6
grp
table(pop(Taiwanmsat2), grp$grp)

dapc.Taiwan <- dapc(Taiwanmsat2, grp$grp)
#keep 15 PCs and all dicriminant functions

dapc.Taiwan
temp<-optim.a.score(dapc.Taiwan) #best # is 5
scatter(dapc.Taiwan, cstar=0)
dapc.Taiwan$posterior

#col=c("blue", "black", "forestgreen", "red", "gray", "magenta", "gold", "skyblue", "purple")

compoplot(dapc.Taiwan)
```


```{r map}
world <- map_data("world")
Taiwan<-subset(world, region %in% c("China", "Taiwan"))

Taiwan<-ggplot() + 
  geom_polygon(data = Taiwan, aes(x=long, y = lat, group = group), fill = NA, color = "black") +
  coord_map(xlim = c(115,124),ylim = c(21,27.5))#, 1.3)

#Taiwan<-ggplot() + 
#  geom_polygon(data = Taiwan, aes(x=long, y = lat, group = group), fill = NA, color = "black") + 
#  coord_fixed(1.3)
GPS_Taiwan<-read.table("C:/Users/Sandra/OneDrive/kudzu/Taiwan/least_missing_GPS.txt", header = T)
GPS_Taiwan
length(GPS_Taiwan[,1])

labs <- data.frame(
  long=GPS_Taiwan[,3],
  lat=GPS_Taiwan[,2],
  names=GPS_Taiwan[,1],
  stringsAsFactors = FALSE
); labs

data <-as.data.frame(unclass(table(pop(Taiwanmsat2), grp$grp))); data
write.csv(as.matrix(Fst), "C:/Users/Sandra/OneDrive/kudzu/Taiwan/Tw.Ch_Gst_msat.csv")
colnames(data)<- c("Group 1","Group 2","Group 3","Group 4") #,"Group 5")#,"Group 6")#,"Group 7")

Taiwan_data<-cbind(GPS_Taiwan[,2:3], data); Taiwan_data

#put sampling points on the map
Taiwan + 
  geom_point(data = labs, aes(x = long, y = lat), color = "black", size =2) 

Taiwan +
  geom_scatterpie(data=Taiwan_data, aes(x=Long, y=Lat, r=.08),sorted_by_radius = F, cols= c("Group 1","Group 2","Group 3","Group 4")) + #,"Group 5")) + #,"Group 6")) + #"Group 5","Group 6","Group 7")) +
  scale_fill_manual(values=rainbow(4))


ggsave("C:/Users/Sandra/OneDrive/kudzu/Taiwan/least_missing_GPS_jan28_2021_4grp.pdf", device="pdf", width = 20, height = 20, units ="cm", dpi=600)

```