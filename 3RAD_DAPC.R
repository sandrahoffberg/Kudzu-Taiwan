library(poppr)
devtools::install_github("dwinter/mmod", ref="bs_fix")
library(mmod)
library(adegenet)
library(vcfR)
library(dartR)

vcf = read.vcfR("C:/Users/Sandra/OneDrive/kudzu/3RAD_analysis/kudzu91_noHC_min20_100Loci.vcf")

my_genlight<-vcfR2genlight(vcf)
my_genlight

#number of loci
nLoc(my_genlight)
#number of individials
nInd(my_genlight)


pop(my_genlight) <- as.factor(c(
  "AL12",
  "AL12",
  "AL12",
  "AL12",
  "AL12",
  "AL12",
  "AL12",
  "AL12",
  "AL12",
  "AL12",
  "AL5",
  "AL5",
  "AL5",
  "AL5",
  "AL5",
  "AL5",
  "AL5",
  "AL5",
  "AL5",
  "AL5",
  "AR2",
  "AR2",
  "AR2",
  "AR2",
  "AR2",
  "AR2",
  "FL3",
  "FL3",
  "FL3",
  "FL3",
  "FL3",
  "GA3",
  "GA3",
  "GA3",
  "GA3",
  "GA3",
  "GA3",
  "GA3",
  "GA3",
  "GA34",
  "GA34",
  "GA34",
  "GA34",
  "GA34",
  "GA34",
  "GA34",
  "GA34",
  "GA34",
  "GA34",
  "GA36",
  "GA36",
  "GA36",
  "GA36",
  "GA36",
  "GA36",
  "GA96",
  "GA96",
  "GA96",
  "GA96",
  "GA96",
  "GA96",
  "GA96",
  "GA96",
  "GA96",
  "GA96",
  "GA96",
  "GA96",
  "KAN1",
  "KAN1",
  "KAN1",
  "KAN1",
  "KAN1",
  "KAN4",
  "KAN4",
  "KAN4",
  "KAN5",
  "KAN5",
  "KAN5",
  "KAN5",
  "KAN5",
  "KB1",
  "KB1",
  "KB1",
  "KB1",
  "KB1",
  "KFU7",
  "KFU7",
  "KFU7",
  "KFU7",
  "KFU7",
  "KFU7",
  "KFU7",
  "KFU7",
  "KFU7",
  "KFU7",
  "KGD3",
  "KGD3",
  "KGD3",
  "KGD3",
  "KGD3",
  "KGD3",
  "KGD3",
  "KGD3",
  "KGD3",
  "KGD3",
  "NA",
  "KGD4",
  "KGD4",
  "NA",
  "NA",
  "KGD4",
  "KGD4",
  "NA",
  "KGD4",
  "NA",
  "NA",
  "KGD5",
  "KGD5",
  "KGU1",
  "KGU1",
  "KGU1",
  "KGU1",
  "KGU1",
  "KGU1",
  "KGU1",
  "KGU1",
  "KGU1",
  "KGU1",
  "KGU1",
  "KHN10",
  "KHN10",
  "KHN10",
  "KHN10",
  "KHN10",
  "KHN10",
  "KHN10",
  "KHN10",
  "KHN10",
  "KHN10",
  "KHN11",
  "KHN11",
  "NA",
  "KHN11",
  "KHN5",
  "KHN5",
  "KHN5",
  "KHN5",
  "KHN5",
  "KHN5",
  "KHN5",
  "KHN5",
  "KHN5",
  "KHN5",
  "KHN7",
  "KHN7",
  "KHN7",
  "KHN7",
  "KHN7",
  "KHN7",
  "KHN7",
  "KHN7",
  "KHN7",
  "KHN7",
  "KHN9",
  "KHN9",
  "KHN9",
  "KHN9",
  "KHN9",
  "KHN9",
  "KHN9",
  "KHN9",
  "KHN9",
  "KHN9",
  "KHU3",
  "KHU2",
  "KHU2",
  "KHU2",
  "KHU2",
  "KHU2",
  "KHU2",
  "KHU2",
  "KHU2",
  "KHU2",
  "KHU3",
  "KHU3",
  "KHU3",
  "KHU3",
  "KHU3",
  "KHU3",
  "KHU3",
  "KHU3",
  "KHU3",
  "KJI2",
  "KJI2",
  "KJI2",
  "KJI2",
  "KJI2",
  "KJI2",
  "KJI2",
  "KJI2",
  "KJI2",
  "KJI2",
  "KJI5",
  "KJI5",
  "KJI5",
  "KJI5",
  "KJI5",
  "KJIX2",
  "KJIX2",
  "KJIX2",
  "KJIX2",
  "KJIX2",
  "KJIX2",
  "KJIX2",
  "KJIX2",
  "KJIX2",
  "KJIX2",
  "KJP11",
  "KJP11",
  "KJP11",
  "KJP11",
  "KJP11",
  "KJP11",
  "KJP11",
  "KJP11",
  "KJP11",
  "KJP11",
  "KJP14",
  "KJP14",
  "KJP14",
  "KJP14",
  "KJP14",
  "KJP14",
  "KJP14",
  "KJP14",
  "KJP14",
  "KJP14",
  "KJP17",
  "KJP17",
  "KJP17",
  "KJP17",
  "KJP17",
  "KJP18",
  "KJP18",
  "KJP18",
  "KJP18",
  "KJP18",
  "KJP18",
  "KJP18",
  "KJP18",
  "KJP18",
  "KJP18",
  "KJP2",
  "KJP2",
  "KJP2",
  "KJP2",
  "KJP2",
  "KJP2",
  "KJP2",
  "KJP2",
  "KJP2",
  "KJP2",
  "KJP22",
  "KJP22",
  "KJP22",
  "KJP22",
  "KJP22",
  "KJP22",
  "KJP22",
  "KJP22",
  "KJP22",
  "KJP22",
  "KJP23",
  "KJP23",
  "KJP23",
  "KJP23",
  "KJP23",
  "KJP23",
  "KJP23",
  "KJP23",
  "KJP23",
  "KJP23",
  "KJP25",
  "KJP25",
  "KJP25",
  "KJP25",
  "KJP25",
  "KJP25",
  "KJP25",
  "KJP25",
  "KJP25",
  "KJP25",
  "KJP27",
  "KJP27",
  "KJP27",
  "KJP27",
  "KJP27",
  "KJP27",
  "KJP27",
  "KJP27",
  "KJP27",
  "KJP3",
  "KJP3",
  "KJP3",
  "KJP3",
  "KJP3",
  "KJP3",
  "KJP3",
  "KJP3",
  "KJP3",
  "KJP3",
  "KJP5",
  "KJP5",
  "KJP5",
  "KJP5",
  "KJP5",
  "KJP5",
  "KJP5",
  "KJP5",
  "KJP5",
  "KJP5",
  "KJP8",
  "KJP8",
  "KJP8",
  "KJP8",
  "KJP8",
  "KJP8",
  "KJP8",
  "KJP8",
  "KJP8",
  "KJP8",
  "KKO5",
  "KKO19",
  "KKO19",
  "KKO19",
  "KKO19",
  "KKO19",
  "KKO19",
  "KKO19",
  "KKO19",
  "KKO19",
  "KKO19",
  "KKO19",
  "KKO22",
  "KKO22",
  "KKO22",
  "KKO22",
  "KKO22",
  "KKO22",
  "KKO22",
  "KKO22",
  "KKO22",
  "KKO22",
  "KKO22",
  "KKO3",
  "KKO3",
  "KKO3",
  "KKO3",
  "KKO3",
  "KKO3",
  "KKO3",
  "KKO3",
  "KKO3",
  "KKO3",
  "KKO5",
  "KKO5",
  "KKO5",
  "KKO5",
  "KKO5",
  "KKO5",
  "KKO5",
  "KKO5",
  "KKO5",
  "KKO6",
  "KKO6",
  "KKO6",
  "KKO6",
  "KKO6",
  "KKO6",
  "KKO6",
  "KKO6",
  "KKO6",
  "KKO6",
  "KLN1",
  "KLN1",
  "KLN1",
  "KLN2",
  "KLN2",
  "KLN2",
  "KLN2",
  "KLN2",
  "KSA1",
  "KSA1",
  "KSA1",
  "KSA1",
  "KSA1",
  "KSA1",
  "KSA1",
  "KSA1",
  "KSA1",
  "KSA1",
  "KSA2",
  "KSA2",
  "KSA2",
  "KSA2",
  "KSA2",
  "KSA2",
  "KSA2",
  "KSA2",
  "KSA2",
  "KSA2",
  "KSA3",
  "KSA3",
  "KSA3",
  "KSA3",
  "KSA3",
  "KSA3",
  "KSA3",
  "KSA3",
  "KSA3",
  "KSA3",
  "KSC1",
  "KSC1",
  "KSC1",
  "KSC1",
  "KSC1",
  "KSC1",
  "KSC1",
  "KSC1",
  "KSC1",
  "KSC1",
  "KSC4",
  "KSC4",
  "KSC4",
  "KSC4",
  "KSC4",
  "KSC4",
  "KSC4",
  "KSC4",
  "KSC4",
  "KSC4",
  "KSD1",
  "KSD1",
  "KSD1",
  "KSD1",
  "KSD1",
  "KSD1",
  "KSD1",
  "KSD1",
  "KSD1",
  "KSD1",
  "KSD2",
  "KSD2",
  "KSD2",
  "KSD2",
  "KTW1",
  "KTW1",
  "KTW1",
  "KTW1",
  "KTW1",
  "KTW1",
  "KTW1",
  "KTW1",
  "KTW1",
  "KTW1",
  "KTW4",
  "KTW4",
  "KTW4",
  "KTW4",
  "KTW4",
  "KTW4",
  "KTW4",
  "KTW4",
  "KTW4",
  "KTW4",
  "KTW6",
  "KTW6",
  "KTW6",
  "KTW6",
  "KTW6",
  "KTW6",
  "KTW6",
  "KTW6",
  "KTW6",
  "KTW6",
  "KY11",
  "KY11",
  "KY11",
  "KY11",
  "KY11",
  "KY11",
  "KY11",
  "KY11",
  "KY11",
  "KY11",
  "KY7",
  "KY7",
  "KY7",
  "KY7",
  "KY7",
  "KY8",
  "KY8",
  "KY8",
  "KY8",
  "KY8",
  "KY8",
  "KY8",
  "KY8",
  "KY8",
  "KY8",
  "KYU2",
  "KYU2",
  "KYU2",
  "KYU2",
  "KYU2",
  "KYU2",
  "KYU2",
  "KYU2",
  "KYU2",
  "KYU2",
  "KYU7",
  "KYU7",
  "KYU7",
  "KYU7",
  "KYU7",
  "KYU7",
  "KYU7",
  "KYU7",
  "KYU7",
  "KYU7",
  "KZH10",
  "KZH10",
  "KZH10",
  "KZH10",
  "KZH10",
  "KZH10",
  "KZH10",
  "KZH10",
  "KZH10",
  "KZH10",
  "KZH8",
  "KZH8",
  "KZH8",
  "KZH8",
  "KZH8",
  "KZH8",
  "KZH8",
  "KZH8",
  "KZH8",
  "KZH8",
  "MD1",
  "MD1",
  "MD1",
  "MD1",
  "MD1",
  "MD1",
  "MD1",
  "MD1",
  "MD1",
  "MD1",
  "MO1",
  "MO1",
  "MO1",
  "MS4",
  "MS4",
  "MS4",
  "MS4",
  "MS4",
  "MS4",
  "MS4",
  "MS4",
  "MS4",
  "MS4",
  "MS6",
  "MS6",
  "MS6",
  "MS6",
  "MS6",
  "MS6",
  "MS6",
  "MS6",
  "NC13",
  "NC13",
  "NC13",
  "NC13",
  "NC13",
  "NC13",
  "NC13",
  "NC13",
  "NC13",
  "NC20",
  "NC20",
  "NC20",
  "NC20",
  "NC20",
  "NC20",
  "NC20",
  "NC20",
  "NC20",
  "NC20",
  "NC20",
  "NC21",
  "NC21",
  "NC21",
  "NC21",
  "NC21",
  "NC3",
  "NC3",
  "NC3",
  "NC3",
  "NC3",
  "NC3",
  "NC3",
  "NC3",
  "NC3",
  "NC3",
  "NC6",
  "NC6",
  "NC6",
  "NC6",
  "NC6",
  "NC6",
  "NC6",
  "NC6",
  "NC6",
  "NE1",
  "NE1",
  "NE1",
  "NE1",
  "NE1",
  "NI1",
  "NI2",
  "NY4",
  "NY4",
  "NY4",
  "NY4",
  "NY4",
  "OK1",
  "OK1",
  "OK1",
  "OK1",
  "OK1",
  "SC4",
  "SC4",
  "SC4",
  "SC4",
  "SC4",
  "SC4",
  "SC4",
  "SC4",
  "SC4",
  "SC4",
  "SC4",
  "TN6",
  "TN6",
  "TN6",
  "TN6",
  "TN6",
  "TN6",
  "TN6",
  "TN8",
  "TN8",
  "TN8",
  "TN8",
  "TN8",
  "TN8",
  "TN8",
  "TN8",
  "TN8",
  "TN8",
  "TN8",
  "TX4",
  "TX4",
  "TX4",
  "TX4",
  "TX4",
  "TX4",
  "TX4",
  "TX4",
  "TX4",
  "TX4",
  "TX4",
  "TX4",
  "TX5",
  "TX5",
  "TX5",
  "TX5",
  "TX5",
  "TX5",
  "TX5",
  "TX5",
  "TX5",
  "TX5",
  "WV1",
  "WV1",
  "WV1",
  "WV1",
  "WV1",
  "WV1",
  "WV1",
  "WV1",
  "WV1"
  ))

popNames(my_genlight)


#number of populations
nPop(my_genlight)


TW_genlight = popsub(my_genlight, sublist = c("KTW1","KTW4","KTW6","KFU7","KGD3","KGD4","KGD5", "KGU1",
                                "KHN5","KHN9","KHN7","KHN10","KHN11"), drop=TRUE)


#number of individials
nInd(TW_genlight)
TW_genlight@ind.names

gl2<-gl.filter.monomorphs(TW_genlight, verbose=5)
#gl3<-gl.filter.callrate(gl2, method="ind", threshold=0.95, mono.rm=TRUE, plot=TRUE, verbose=5)

gl2

my_genind <- gl2gi(gl2) # not enough memory 

my_genind

save(my_genind, file= "C:/Users/Sandra/OneDrive/kudzu/Taiwan/TW_genind.rdata")

gi2<- missingno(my_genind, type="loci", cutoff = .4)
gi2
gi3<- missingno(my_genind, type="loci", cutoff = .1)
gi3

gi4<-missingno(gi3, type="geno", cutoff=.25)
gi4
#Found 33904 missing values.
#1 genotype contained missing values greater than 25%
#Removing 1 genotype: , KHN10-7A


psex(gi4)

# Now we can collect information of the thresholds. We can set threshold = 1 because we know that this will capture
# the maximum possible distance:
invisible(thresholds <- mlg.filter(gi4, distance = bitwise.dist, stats = "THRESHOLDS", 
                                   threshold = 1, algorithm = "average_neighbor", missing="mean"))

# We can use these thresholds to find an appropriate cutoff
pcut <- cutoff_predictor(thresholds)
pcut
{hist(thresholds, nclass=200)
  abline(v=pcut, col="red")
  abline(v=.035, col="blue")}


```{r DAPC + map}

library(ggplot2)
library(dplyr)
library(stringr)
library(maps)
library(mapdata)
library(mapplots)

library(scatterpie)



#CHoose number of clusters 

set.seed(021921)


###############K=2###############


grp<-find.clusters(gi4) #keep all PCs
150
2

#Keep all PCS and then keep # clusters (elbow)

grp
table(pop(gi4), grp$grp)
table.value(table(pop(gi4), grp$grp), col.lab=paste("inf", 1:50),
            row.lab=paste("ori", 1:45))


#Run DAPC

dapc.Taiwan <- dapc(gi4, grp$grp)
30
15
#Do not keep all PCs. Can keep all in barplot if small #
scatter(dapc.Taiwan, scree.da=FALSE, bg="white", pch=20, cell=0, cstar=0, col=rainbow(6), solid=.4,
        cex=3,clab=0, leg=TRUE, txt.leg=paste("Cluster",1:6))

something <- optim.a.score(dapc.Taiwan)
# 7 and up is good





system.time(pramx <- xvalDapc(tab(gi4, NA.method = "mean"), #pop(Taiwanmsats),
                              grp=grp$grp,
                              n.da=5,
                              #n.pca = 6:14,
                              n.rep = 100,
                              parallel = "multicore", ncpus = 4L))
pramx




######## R map with clones ######
world <- map_data("world")
Taiwan<-subset(world, region %in% c("China", "Taiwan"))


Taiwan<-ggplot() + 
  geom_polygon(data = Taiwan, aes(x=long, y = lat, group = group), fill = NA, color = "black") + 
  coord_map(xlim = c(107,123),ylim = c(20,30)) # coord_fixed(1.3)


GPS_Taiwan<-read.table("C:/Users/Sandra/OneDrive/kudzu/Taiwan/3RAD_GPS_bypop.txt", header = T)
GPS_Taiwan

length(GPS_Taiwan[,1])

labs <- data.frame(
  long=GPS_Taiwan[,3],
  lat=GPS_Taiwan[,2],
  names=GPS_Taiwan[,1],
  stringsAsFactors = FALSE
); labs

data <-as.data.frame(unclass(table(pop(gi4), grp$grp))); data
colnames(data)<- c("Group 1","Group 2")


Taiwan_data<-cbind(GPS_Taiwan[,2:3], data); Taiwan_data


Taiwan +
  geom_scatterpie(data=Taiwan_data, aes(x=Long, y=Lat, r=.2),sorted_by_radius = F, cols= colnames(data)) +
  scale_fill_manual(values=c("red", "blue")) #rainbow(8))


ggsave("C:/Users/Sandra/OneDrive/kudzu/Taiwan/Taiwan_3RAD_map_K2.svg", device="svg", dpi=1200,
       height=6.71, width=10, units="in")



###############K=3###############

#CHoose number of clusters 
#compare between clone corrected and full dataset

grp<-find.clusters(gi4) #keep all PCs
150
3

#Keep all PCS and then keep # clusters (elbow)

grp
table(pop(gi4), grp$grp)
table.value(table(pop(gi4), grp$grp), col.lab=paste("inf", 1:50),
            row.lab=paste("ori", 1:45))


#Run DAPC

dapc.Taiwan <- dapc(gi4, grp$grp)
30
15
#Do not keep all PCs. Can keep all in barplot if small #
scatter(dapc.Taiwan, scree.da=FALSE, bg="white", pch=20, cell=0, cstar=0, col=rainbow(6), solid=.4,
        cex=3,clab=0, leg=TRUE, txt.leg=paste("Cluster",1:6))

something <- optim.a.score(dapc.Taiwan)
# 7 and up is good




system.time(pramx <- xvalDapc(tab(gi4, NA.method = "mean"), #pop(Taiwanmsats),
                              grp=grp$grp,
                              n.da=5,
                              #n.pca = 6:14,
                              n.rep = 100,
                              parallel = "snow", ncpus = 4L))
pramx




######## R map with clones ######
world <- map_data("world")
Taiwan<-subset(world, region %in% c("China", "Taiwan"))


Taiwan<-ggplot() + 
  geom_polygon(data = Taiwan, aes(x=long, y = lat, group = group), fill = NA, color = "black") + 
  coord_map(xlim = c(107,123),ylim = c(20,30)) # coord_fixed(1.3)


GPS_Taiwan<-read.table("C:/Users/Sandra/OneDrive/kudzu/Taiwan/3RAD_GPS_bypop.txt", header = T)
GPS_Taiwan

length(GPS_Taiwan[,1])

labs <- data.frame(
  long=GPS_Taiwan[,3],
  lat=GPS_Taiwan[,2],
  names=GPS_Taiwan[,1],
  stringsAsFactors = FALSE
); labs

data <-as.data.frame(unclass(table(pop(gi4), grp$grp))); data
colnames(data)<- c("Group 1","Group 2","Group 3")


Taiwan_data<-cbind(GPS_Taiwan[,2:3], data); Taiwan_data


Taiwan +
  geom_scatterpie(data=Taiwan_data, aes(x=Long, y=Lat, r=.2),sorted_by_radius = F, cols= colnames(data)) +
  scale_fill_manual(values=c("blue","orange", "red")) #rainbow(8))

ggsave("C:/Users/Sandra/OneDrive/kudzu/Taiwan/Taiwan_3RAD_map_K3.svg", device="svg", dpi=1200,
       height=6.71, width=10, units="in")



###############K=4###############

#CHoose number of clusters 
#compare between clone corrected and full dataset

grp<-find.clusters(gi4) #keep all PCs
150
4

#Keep all PCS and then keep # clusters (elbow)

grp
table(pop(gi4), grp$grp)
table.value(table(pop(gi4), grp$grp), col.lab=paste("inf", 1:50),
            row.lab=paste("ori", 1:45))


#Run DAPC

dapc.Taiwan <- dapc(gi4, grp$grp)
30
15
#Do not keep all PCs. Can keep all in barplot if small #
scatter(dapc.Taiwan, scree.da=FALSE, bg="white", pch=20, cell=0, cstar=0, col=rainbow(6), solid=.4,
        cex=3,clab=0, leg=TRUE, txt.leg=paste("Cluster",1:6))

something <- optim.a.score(dapc.Taiwan)
# 7 and up is good




system.time(pramx <- xvalDapc(tab(gi4, NA.method = "mean"), #pop(Taiwanmsats),
                              grp=grp$grp,
                              n.da=5,
                              #n.pca = 6:14,
                              n.rep = 100,
                              parallel = "snow", ncpus = 4L))
pramx




######## R map with clones ######
world <- map_data("world")
Taiwan<-subset(world, region %in% c("China", "Taiwan"))


Taiwan<-ggplot() + 
  geom_polygon(data = Taiwan, aes(x=long, y = lat, group = group), fill = NA, color = "black") + 
  coord_map(xlim = c(107,123),ylim = c(20,30)) # coord_fixed(1.3)


GPS_Taiwan<-read.table("C:/Users/Sandra/OneDrive/kudzu/Taiwan/3RAD_GPS_bypop.txt", header = T)
GPS_Taiwan

length(GPS_Taiwan[,1])

labs <- data.frame(
  long=GPS_Taiwan[,3],
  lat=GPS_Taiwan[,2],
  names=GPS_Taiwan[,1],
  stringsAsFactors = FALSE
); labs

data <-as.data.frame(unclass(table(pop(gi4), grp$grp))); data
colnames(data)<- c("Group 1","Group 2","Group 3","Group 4")#,"Group 5","Group 6")#,"Group 7",



Taiwan_data<-cbind(GPS_Taiwan[,2:3], data); Taiwan_data


Taiwan +
  geom_scatterpie(data=Taiwan_data, aes(x=Long, y=Lat, r=.2),sorted_by_radius = F, cols= colnames(data)) +
  scale_fill_manual(values=c("orange", "blue", "maroon4", "red")) #rainbow(8))

ggsave("C:/Users/Sandra/OneDrive/kudzu/Taiwan/Taiwan_3RAD_map_K4.svg", device="svg", dpi=1200,
       height=6.71, width=10, units="in")



###############K=5###############

#CHoose number of clusters 
#compare between clone corrected and full dataset

grp<-find.clusters(gi4) #keep all PCs
150
5

#Keep all PCS and then keep # clusters (elbow)

grp
table(pop(gi4), grp$grp)
table.value(table(pop(gi4), grp$grp), col.lab=paste("inf", 1:50),
            row.lab=paste("ori", 1:45))


#Run DAPC

dapc.Taiwan <- dapc(gi4, grp$grp)
30
15
#Do not keep all PCs. Can keep all in barplot if small #
scatter(dapc.Taiwan, scree.da=FALSE, bg="white", pch=20, cell=0, cstar=0, col=rainbow(6), solid=.4,
        cex=3,clab=0, leg=TRUE, txt.leg=paste("Cluster",1:6))

something <- optim.a.score(dapc.Taiwan)
# 7 and up is good




system.time(pramx <- xvalDapc(tab(gi4, NA.method = "mean"), #pop(Taiwanmsats),
                              grp=grp$grp,
                              n.da=5,
                              #n.pca = 6:14,
                              n.rep = 100,
                              parallel = "multicore", ncpus = 4L))
pramx




######## R map with clones ######
world <- map_data("world")
Taiwan<-subset(world, region %in% c("China", "Taiwan"))


Taiwan<-ggplot() + 
  geom_polygon(data = Taiwan, aes(x=long, y = lat, group = group), fill = NA, color = "black") + 
  coord_map(xlim = c(107,123),ylim = c(20,30)) # coord_fixed(1.3)


GPS_Taiwan<-read.table("C:/Users/Sandra/OneDrive/kudzu/Taiwan/3RAD_GPS_bypop.txt", header = T)
GPS_Taiwan

length(GPS_Taiwan[,1])

labs <- data.frame(
  long=GPS_Taiwan[,3],
  lat=GPS_Taiwan[,2],
  names=GPS_Taiwan[,1],
  stringsAsFactors = FALSE
); labs

data <-as.data.frame(unclass(table(pop(gi4), grp$grp))); data
colnames(data)<- c("Group 1","Group 2","Group 3","Group 4","Group 5")


Taiwan_data<-cbind(GPS_Taiwan[,2:3], data); Taiwan_data


Taiwan +
  geom_scatterpie(data=Taiwan_data, aes(x=Long, y=Lat, r=.2),sorted_by_radius = F, cols= colnames(data)) +
  scale_fill_manual(values=c("blue", "hotpink", "red", "maroon4", "orange")) #rainbow(8))

ggsave("C:/Users/Sandra/OneDrive/kudzu/Taiwan/Taiwan_3RAD_map_K5.svg", device="svg", dpi=1200,
       height=6.71, width=10, units="in")






###############K=6###############

#CHoose number of clusters 
#compare between clone corrected and full dataset

grp<-find.clusters(gi4) #keep all PCs
150
6

#Keep all PCS and then keep # clusters (elbow)

grp
table(pop(gi4), grp$grp)
table.value(table(pop(gi4), grp$grp), col.lab=paste("inf", 1:50),
            row.lab=paste("ori", 1:45))


#Run DAPC

dapc.Taiwan <- dapc(gi4, grp$grp)
30
15
#Do not keep all PCs. Can keep all in barplot if small #
scatter(dapc.Taiwan, scree.da=FALSE, bg="white", pch=20, cell=0, cstar=0, col=rainbow(6), solid=.4,
        cex=3,clab=0, leg=TRUE, txt.leg=paste("Cluster",1:6))

something <- optim.a.score(dapc.Taiwan)
# 7 and up is good





system.time(pramx <- xvalDapc(tab(gi4, NA.method = "mean"), #pop(Taiwanmsats),
                              grp=grp$grp,
                              n.da=5,
                              #n.pca = 6:14,
                              n.rep = 100,
                              parallel = "multicore", ncpus = 4L))
pramx




##### R map with clones #####
world <- map_data("world")
Taiwan<-subset(world, region %in% c("China", "Taiwan"))


Taiwan<-ggplot() + 
  geom_polygon(data = Taiwan, aes(x=long, y = lat, group = group), fill = NA, color = "black") + 
  coord_map(xlim = c(107,123),ylim = c(20,30)) # coord_fixed(1.3)


GPS_Taiwan<-read.table("C:/Users/Sandra/OneDrive/kudzu/Taiwan/3RAD_GPS_bypop.txt", header = T)
GPS_Taiwan

length(GPS_Taiwan[,1])

labs <- data.frame(
  long=GPS_Taiwan[,3],
  lat=GPS_Taiwan[,2],
  names=GPS_Taiwan[,1],
  stringsAsFactors = FALSE
); labs

data <-as.data.frame(unclass(table(pop(gi4), grp$grp))); data
colnames(data)<- c("Group 1","Group 2","Group 3","Group 4","Group 5","Group 6")


Taiwan_data<-cbind(GPS_Taiwan[,2:3], data); Taiwan_data


Taiwan +
  geom_scatterpie(data=Taiwan_data, aes(x=Long, y=Lat, r=.2),sorted_by_radius = F, cols= colnames(data)) +
  scale_fill_manual(values=c("sienna", "maroon4", "hotpink", "orange","blue","red")) #rainbow(8))

ggsave("C:/Users/Sandra/OneDrive/kudzu/Taiwan/Taiwan_3RAD_map_K6.svg", device="svg", dpi=1200,
       height=6.71, width=10, units="in")
