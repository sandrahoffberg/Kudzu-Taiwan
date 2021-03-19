devtools::install_github("dkahle/ggmap")
#https://eriqande.github.io/rep-res-web/lectures/making-maps-with-R.html
library(ape)
library(ggplot2)
library(devtools)
library(dplyr)
library(stringr)
library(maps)
library(mapdata)
library(mapplots)
library(scatterpie)
library(svglite)



########## China ######################
input <- "D:/kudzu2019/China/cpDNA_China_good2.fasta"
d <- ape::read.dna(input, format='fasta')
e <- dist.dna(d)
h <- pegas::haplotype(d)
(net <- pegas::haploNet(h))


ind.hap<-with(
  stack(setNames(attr(h, "index"), rownames(h))),
  table(hap=ind, pop=rownames(d)[values])
)

length(ind.hap[,1])

#There are 82 populations
colors82<-rainbow(82, s = 1, v = 1, start = 0, end = max(1, 82 - 1)/82,
                  alpha = 1)#, rev = FALSE)

plot(net, size=attr(net, "freq")/50, scale.ratio = .5, pie=ind.hap, 
     lwd=2, cex=0.8, legend=F, fast=T, show.mutation = 1, threshold = 0, bg=colors82)
legend(-29,15, colnames(ind.hap), col=colors82, pch=19, cex=0.4)






#See the geographic distribution by plotting these on a map

world <- map_data("world")
China<-subset(world, region %in% c("China"))


China<-ggplot() + 
  geom_polygon(data = China, aes(x=long, y = lat, group = group), fill = NA, color = "black") + 
  coord_fixed(1.3)
GPS_China<-read.table("D:/kudzu2019/China/China_cpDNA_GPS_bypop.txt", header = T)
GPS_China

length(GPS_China[,1])

labs <- data.frame(
  long=GPS_China[,3],
  lat=GPS_China[,2],
  names=GPS_China[,1],
  stringsAsFactors = FALSE
)  



labs


#If you want to keep all the data the same but just zoom in, you can use the xlim and ylim arguments to coord_cartesian(). Though, to keep the aspect ratio correct we must use coord_fixed() instead of coord_cartesian().
#This chops stuff off but doesn’t discard it from the data set:

#  eb2 + coord_fixed(xlim = c(-123, -121.0),  ylim = c(36, 38), ratio = 1.3)


colors20<-rainbow(20)#, s = 1, v = 1, start = 0, end = max(1, 20 - 1)/20, alpha = 1)#, rev = FALSE)



ind2<- unclass(t(ind.hap))
colnames(ind2)<- c("Haplotype I","Haplotype II","Haplotype III","Haplotype IV","Haplotype V",
                   "Haplotype VI","Haplotype VII","Haplotype VIII","Haplotype IX","Haplotype X",
                   "Haplotype XI","Haplotype XII","Haplotype XIII","Haplotype XIV")
ind2


China_data<-cbind(GPS_China[,2:3], ind2)
China_data


#put sampling points on the map
China + 
  geom_point(data = labs, aes(x = long, y = lat), color = "black", size =2) 



China +
  geom_scatterpie(data=China_data, aes(x=Long, y=Lat, r= .7),sorted_by_radius = F, cols= c("Haplotype I","Haplotype II","Haplotype III","Haplotype IV","Haplotype V",
                                                                               "Haplotype VI","Haplotype VII","Haplotype VIII","Haplotype IX","Haplotype X",
                                                                               "Haplotype XI","Haplotype XII","Haplotype XIII","Haplotype XIV")) +
  scale_fill_manual(values=c("black","cyan","white","gray","red","yellow","darkorange","green","darkgreen","skyblue2","purple","blue","hotpink","gold"))



##Now separate the 2 groups in China

input <- "D:/kudzu2019/cpDNA_China_good3.fasta"
d <- ape::read.dna(input, format='fasta')
e <- dist.dna(d)
h <- pegas::haplotype(d)
(net <- pegas::haploNet(h))


ind.hap<-with(
  stack(setNames(attr(h, "index"), rownames(h))),
  table(hap=ind, pop=rownames(d)[values])
)

col2<-c("purple4", "maroon")
plot(net, size=attr(net, "freq")/50, scale.ratio = .5, pie=ind.hap, 
     lwd=2, cex=0.8, legend=F, fast=T, show.mutation = 1, threshold = 0, bg=col2)
legend(-28,0, colnames(ind.hap), col=col2, pch=19, cex=1.3)



#You wouldn't really want to plot this on a map because you are using PCA data to separate groups, 
#not cp data. 




########## Taiwan ######################

#only 1 group, no need to  do it by pop

input <- "C:/Users/Sandra/OneDrive/kudzu/Taiwan/cpDNA_Taiwan_by_pop_noKYU.fasta"
d <- ape::read.dna(input, format='fasta')
e <- dist.dna(d)
h <- pegas::haplotype(d)
(net <- pegas::haploNet(h))


ind.hap<-with(
  stack(setNames(attr(h, "index"), rownames(h))),
  table(hap=ind, pop=rownames(d)[values])
)

sum(ind.hap[1,])
sum(ind.hap[2,])
sum(ind.hap[3,])
sum(ind.hap[4,])
sum(ind.hap[5,])
sum(ind.hap[6,])
sum(ind.hap[7,])
sum(ind.hap[8,])
sum(ind.hap[9,])
sum(ind.hap[1,])/sum(ind.hap)

write.csv(ind.hap, "C:/Users/Sandra/OneDrive/kudzu/Taiwan/Taiwan_cpDNA_haplotypes_bypop.csv")


svg("C:/Users/Sandra/OneDrive/kudzu/Taiwan/Taiwan_cpDNA_haplotypes.svg")

plot(net, size=attr(net, "freq")/100, scale.ratio = 1, pie=ind.hap, 
     lwd=5, cex=0.8, legend=F, fast=T, show.mutation = 1, threshold = 0, #col = c("red", "blue", "black"),
     bg=c("gray")
#,"gray", "red", "blue"))#,"white","white","white","white","white","white","white","white","white","white","white"))
#legend(1,1, colnames(ind.hap), col=c("black","gray", "red", "blue", "forestgreen","gold", "dodgerblue", "deeppink"), pch=19, cex=1)

dev.off()




#See the geographic distribution by plotting these on a map

world <- map_data("world")
TW<-subset(world, region %in% c("China", "Taiwan"))

TW<-ggplot() + 
  geom_polygon(data = TW, aes(x=long, y = lat, group = group), fill = NA, color = "black") +
  coord_map(xlim = c(107,123),ylim = c(20,30))#, 1.3)

TW
GPS_TW<-read.table("C:/Users/Sandra/OneDrive/kudzu/Taiwan/GPS_by_pop.txt", header = T)
GPS_TW

length(GPS_TW[,1])

labs <- data.frame(
  long=GPS_TW[,3],
  lat=GPS_TW[,2],
  names=GPS_TW[,1],
  stringsAsFactors = FALSE
)  



labs


#If you want to keep all the data the same but just zoom in, you can use the xlim and ylim arguments to coord_cartesian(). Though, to keep the aspect ratio correct we must use coord_fixed() instead of coord_cartesian().
#This chops stuff off but doesn’t discard it from the data set:

#  eb2 + coord_fixed(xlim = c(-123, -121.0),  ylim = c(36, 38), ratio = 1.3)


colors3<-rainbow(9)#, s = 1, v = 1, start = 0, end = max(1, 20 - 1)/20, alpha = 1)#, rev = FALSE)



ind2<- unclass(t(ind.hap))
colnames(ind2)<- c("Haplotype I","Haplotype II","Haplotype III",
                   "Haplotype IV","Haplotype V","Haplotype VI",
                   "Haplotype VII","Haplotype VIII","Haplotype IX")
ind2


TW_data<-cbind(GPS_TW[,2:3], ind2)
TW_data


#put sampling points on the map
TW + 
  geom_point(data = labs, aes(x = long, y = lat), color = "black", size =2) 



TW +
  geom_scatterpie(data=TW_data, aes(x=Long, y=Lat, r= .2),sorted_by_radius = F, cols= c("Haplotype I","Haplotype II","Haplotype III",
                                                                                        "Haplotype IV","Haplotype V","Haplotype VI",
                                                                                        "Haplotype VII","Haplotype VIII","Haplotype IX")) +
  scale_fill_manual(values=c("blue","red","orange","yellow", "olivedrab3","mediumspringgreen","deepskyblue","darkviolet","deeppink"))


ggsave("C:/Users/Sandra/OneDrive/kudzu/Taiwan/Taiwan_cpDNA_map_noKYU.svg", device="svg", dpi=600)



########## Korea ######################


input <- "D:/kudzu2019/Korea/cpDNA_Korea_good2.fasta"
d <- ape::read.dna(input, format='fasta')
e <- dist.dna(d)
h <- pegas::haplotype(d)
(net <- pegas::haploNet(h))


ind.hap<-with(
  stack(setNames(attr(h, "index"), rownames(h))),
  table(hap=ind, pop=rownames(d)[values])
)

#There are 21 populations
colors21<-rainbow(21, s = 1, v = 1, start = 0, end = max(1, 21 - 1)/21,
                  alpha = 1)#, rev = FALSE)


plot(net, size=attr(net, "freq")/50, scale.ratio = .5, pie=ind.hap, 
     lwd=2, cex=0.8, legend=F, fast=T, show.mutation = 1, threshold = 0,bg=colors21)
legend(-1.5,1, colnames(ind.hap), col=colors21, pch=19, cex=1)




#See the geographic distribution by plotting these on a map

world <- map_data("world")
Korea<-subset(world, region %in% c("South Korea"))
Korea

Korea<-ggplot() + 
  geom_polygon(data = Korea, aes(x=long, y = lat, group = group), fill = NA, color = "black") + 
  coord_fixed(1.3)

GPS_Korea<-read.table("D:/kudzu2019/Korea/Korea_cpDNA_GPS_bypop.txt", header = T)
GPS_Korea

length(GPS_Korea[,1])

labs <- data.frame(
  long=GPS_Korea[,3],
  lat=GPS_Korea[,2],
  names=GPS_Korea[,1],
  stringsAsFactors = FALSE
);labs


#If you want to keep all the data the same but just zoom in, you can use the xlim and ylim arguments to coord_cartesian(). Though, to keep the aspect ratio correct we must use coord_fixed() instead of coord_cartesian().
#This chops stuff off but doesn’t discard it from the data set:

#  eb2 + coord_fixed(xlim = c(-123, -121.0),  ylim = c(36, 38), ratio = 1.3)

length(ind.hap[,1])
#we need as many colors as there are haplotypes
colors6<-rainbow(6)#, s = 1, v = 1, start = 0, end = max(1, 20 - 1)/20, alpha = 1)#, rev = FALSE)



ind2<- unclass(t(ind.hap))
colnames(ind2)<- c("Haplotype I","Haplotype II","Haplotype III","Haplotype IV","Haplotype V",
                   "Haplotype VI")
ind2


Korea_data<-cbind(GPS_Korea[,2:3], ind2)
Korea_data


#put sampling points on the map
Korea + 
  geom_point(data = labs, aes(x = long, y = lat), color = "black", size =2) 



Korea +
  geom_scatterpie(data=Korea_data, aes(x=Long, y=Lat, r= .1),sorted_by_radius = F, cols= colnames(ind2)) +
  scale_fill_manual(values=c("red","orange","yellow","darkgreen","blue","purple")) +
  theme(panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent"))





########## Japan ######################


input <- "D:/kudzu2019/Japan/cpDNA_Japan_good2.fasta"
d <- ape::read.dna(input, format='fasta')
e <- dist.dna(d)
h <- pegas::haplotype(d)
(net <- pegas::haploNet(h))


ind.hap<-with(
  stack(setNames(attr(h, "index"), rownames(h))),
  table(hap=ind, pop=rownames(d)[values])
)

#There are 29 populations
colors29<-rainbow(29, s = 1, v = 1, start = 0, end = max(1, 29 - 1)/29,
                  alpha = 1)#, rev = FALSE)


plot(net, size=attr(net, "freq")/50, scale.ratio = .5, pie=ind.hap, 
     lwd=2, cex=0.8, legend=F, fast=T, show.mutation = 1, threshold = 0, bg=colors29)
legend(-6,5, colnames(ind.hap), col=colors21, pch=19, cex=.85)




#See the geographic distribution by plotting these on a map

world <- map_data("world")
Japan<-subset(world, region %in% c("Japan"))
Japan

Japan<-ggplot() + 
  geom_polygon(data = Japan, aes(x=long, y = lat, group = group), fill = NA, color = "black") + 
  coord_fixed(1.3)

GPS_Japan<-read.table("D:/kudzu2019/Japan/Japan_cpDNA_GPS_bypop.txt", header = T)
GPS_Japan

length(GPS_Japan[,1])

labs <- data.frame(
  long=GPS_Japan[,3],
  lat=GPS_Japan[,2],
  names=GPS_Japan[,1],
  stringsAsFactors = FALSE
);labs


#If you want to keep all the data the same but just zoom in, you can use the xlim and ylim arguments to coord_cartesian(). Though, to keep the aspect ratio correct we must use coord_fixed() instead of coord_cartesian().
#This chops stuff off but doesn’t discard it from the data set:

#  eb2 + coord_fixed(xlim = c(-123, -121.0),  ylim = c(36, 38), ratio = 1.3)

length(ind.hap[,1])
#we need as many colors as there are haplotypes
colors8<-rainbow(8)#, s = 1, v = 1, start = 0, end = max(1, 20 - 1)/20, alpha = 1)#, rev = FALSE)



ind2<- unclass(t(ind.hap))
colnames(ind2)<- c("Haplotype I","Haplotype II","Haplotype III","Haplotype IV","Haplotype V",
                   "Haplotype VI","Haplotype VII","Haplotype VIII")
ind2


Japan_data<-cbind(GPS_Japan[,2:3], ind2)
Japan_data


#put sampling points on the map
Japan + 
  geom_point(data = labs, aes(x = long, y = lat), color = "black", size =2) 



Japan +
  geom_scatterpie(data=Japan_data, aes(x=Long, y=Lat, r= .4),sorted_by_radius = F, cols= colnames(ind2)) +
  scale_fill_manual(values=c("red","orange","yellow","darkgreen","cyan", "blue","purple","magenta")) +
  theme(panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent"))




########## AllAsia ######################


input <- "D:/kudzu2019/cpDNA_AllAsia.fasta"
d <- ape::read.dna(input, format='fasta')
e <- dist.dna(d)
h <- pegas::haplotype(d)
(net <- pegas::haploNet(h))


ind.hap<-with(
  stack(setNames(attr(h, "index"), rownames(h))),
  table(hap=ind, pop=rownames(d)[values])
)


#There are 5 populations
colors4<-c("gold","deeppink","blue","green")


plot(net, size=attr(net, "freq")/40, scale.ratio = .5, pie=ind.hap, 
     lwd=2, cex=0.8, legend=F, fast=T, show.mutation = 1, threshold = 0, bg=colors4)
legend(-15,5, colnames(ind.hap), col=colors4, pch=19, cex=.85)




#this shows the haplotype network. To put it on a map, need populations: 
#AllAsia

input <- "D:/kudzu2019/cpDNA_Asia_by_pop.fasta"
d <- ape::read.dna(input, format='fasta')
e <- dist.dna(d)
h <- pegas::haplotype(d)
(net <- pegas::haploNet(h))


ind.hap<-with(
  stack(setNames(attr(h, "index"), rownames(h))),
  table(hap=ind, pop=rownames(d)[values])
)



ind.hap
ind.hap2<-subset(as.data.frame(as.table(ind.hap)), Freq>0)
ind.hap2
ind.hap3 <- ind.hap2[rep(row.names(ind.hap2), ind.hap2$Freq), 1:2]
ind.hap3
length(ind.hap3[,1])


#write.table(ind.hap2, "D:/kudzu2019/cpDNA_haplotypes_AllAsia_bypop.txt")





###put this on a map





#Plotting maps-package maps with ggplot
# *  The maps package contains a lot of outlines of continents, countries, states, and counties that have been with R for a long time.
# *  The mapdata package contains a few more, higher-resolution outlines.
# *  The maps package comes with a plotting function, but, we will opt to use ggplot2 to plot the maps in the maps package.
# *  Recall that ggplot2 operates on data frames. Therefore we need some way to translate the maps data into a data frame format the ggplot can use.

#Maps in the maps package
#Package maps provides lots of different map outlines and points for cities, etc.
#Some examples: usa, nz, state, world, etc.
#Makin’ data frames from map outlines
#ggplot2 provides the map_data() function.
#Think of it as a function that turns a series of points along an outline into a data frame of those points.
#Syntax: map_data("name") where “name” is a quoted string of the name of a map in the maps or mapdata package
#Here we get a USA map from maps:

world <- map_data("world")

dim(world)
#> [1] 99338    6

head(world)
#>        long      lat group order region subregion
#> 1 -101.4078 29.74224     1     1   main      <NA>
#> 2 -101.3906 29.74224     1     2   main      <NA>
#> 3 -101.3620 29.65056     1     3   main      <NA>
#> 4 -101.3505 29.63911     1     4   main      <NA>
#> 5 -101.3219 29.63338     1     5   main      <NA>
#> 6 -101.3047 29.64484     1     6   main      <NA>


# * order. This just shows in which order ggplot should “connect the dots”
# * region and subregion tell what region or subregion a set of points surrounds.
# * group. This is very important! ggplot2’s functions can take a group argument which controls (amongst other things) whether adjacent points should be connected by lines. If they are in the same group, then they get connected, but if they are in different groups then they don’t.
#Essentially, having to points in different groups means that ggplot “lifts the pen” when going between them.

#Maps in this format can be plotted with the polygon geom. i.e. using geom_polygon().
#geom_polygon() drawn lines between points and “closes them up” (i.e. draws a line from the last point back to the first point)
#You have to map the group aesthetic to the group column
#Of course, x = long and y = lat are the other aesthetics.
#Simple black map
#By default, geom_polygon() draws with no line color, but with a black fill:

world <- map_data("world") # we already did this, but we can do it again
ggplot() + geom_polygon(data = world, aes(x=long, y = lat, group = group)) + 
  coord_fixed(1.3)


#What is this coord_fixed()?
#  This is very important when drawing maps.
#It fixes the relationship between one unit in the y direction and one unit in the x direction.
#Then, even if you change the outer dimensions of the plot (i.e. by changing the window size or the size of the pdf file you are saving it to (in ggsave for example)), the aspect ratio remains unchanged.
#In the above case, I decided that if every y unit was 1.3 times longer than an x unit, then the plot came out looking good.
#A different value might be needed closer to the poles.
#Mess with line and fill colors
#Here is no fill, with a red line. Remember, fixed value of aesthetics go outside the aes function.


Asia<-subset(world, region %in% c("North Korea","China","South Korea", "Japan"))


sort(unique(map_data("world")$region))


Asia<-ggplot() + 
  geom_polygon(data = Asia, aes(x=long, y = lat, group = group), fill = NA, color = "black") + 
  coord_fixed(1.3)



#Adding points to the map

GPS_points<-read.table("D:/kudzu2019/Asia_cpDNA_GPS_bypop.txt", header = T)
GPS_points

length(GPS_points[,1])

labs <- data.frame(
  long=GPS_points[,3],
  lat=GPS_points[,2],
  names=GPS_points[,1],
  stringsAsFactors = FALSE
)  



labs


#If you want to keep all the data the same but just zoom in, you can use the xlim and ylim arguments to coord_cartesian(). Though, to keep the aspect ratio correct we must use coord_fixed() instead of coord_cartesian().
#This chops stuff off but doesn’t discard it from the data set:

#  eb2 + coord_fixed(xlim = c(-123, -121.0),  ylim = c(36, 38), ratio = 1.3)


colors20<-rainbow(20)#, s = 1, v = 1, start = 0, end = max(1, 20 - 1)/20, alpha = 1)#, rev = FALSE)

length(GPS_points[,3])
length(as.numeric(ind.hap))
length(ind.hap3$pop)



ind2<- unclass(t(ind.hap))
colnames(ind2)<- c("Haplotype I","Haplotype II","Haplotype III","Haplotype IV","Haplotype V",
                   "Haplotype VI","Haplotype VII","Haplotype VIII","Haplotype IX","Haplotype X",
                   "Haplotype XI","Haplotype XII","Haplotype XIII","Haplotype XIV","Haplotype XV","Haplotype XVI",
                   "Haplotype XVII","Haplotype XVIII","Haplotype XIX","Haplotype XX")
ind2


yes<-NULL
yes<-cbind(GPS_points[,2:3], ind2)
yes


#put sampling points on the map
Asia + 
  geom_point(data = labs, aes(x = long, y = lat), color = "black", size =2) 



Asia +
  geom_scatterpie(data=yes, aes(x=Y, y=X, r= .7),sorted_by_radius = F, cols= c("Haplotype I","Haplotype II","Haplotype III","Haplotype IV","Haplotype V",
                                                                               "Haplotype VI","Haplotype VII","Haplotype VIII","Haplotype IX","Haplotype X",
                                                                               "Haplotype XI","Haplotype XII","Haplotype XIII","Haplotype XIV","Haplotype XV","Haplotype XVI",
                                                                               "Haplotype XVII","Haplotype XVIII","Haplotype XIX","Haplotype XX")) +
  scale_fill_manual(values=rainbow(21))



# p <- ggplot(Asia, aes(long, lat)) +
#   geom_map(map=Asia, aes(map_id=region), fill=NA, color="black") +
#   coord_quickmap()
# p + geom_scatterpie(aes(x=long, y=lat, group=region, r=radius),
#                     data=d, cols=LETTERS[1:4], color=NA, alpha=.8) +
#   geom_scatterpie_legend(d$radius, x=-160, y=-55)







ggmap
The ggmap package is the most exciting R mapping tool in a long time! You might be able to get better looking maps at some resolutions by using shapefiles and rasters from naturalearthdata.com but ggmap will get you 95% of the way there with only 5% of the work!
  
  Three examples
I am going to run through three examples. Working from the small spatial scale up to a larger spatial scale.
Named “sampling” points on the Sisquoc River from the “Sisquoctober Adventure”
A GPS track from a short bike ride in Wilder Ranch.
Fish sampling locations from the coded wire tag data base.
How ggmap works
ggmap simplifies the process of downloading base maps from Google or Open Street Maps or Stamen Maps to use in the background of your plots.
It also sets the axis scales, etc, in a nice way.
Once you have gotten your maps, you make a call with ggmap() much as you would with ggplot()
Let’s do by example.
Sisquoctober
Here is a small data frame of points from the Sisquoc River.

sisquoc <- read.table("data/sisquoc-points.txt", sep = "\t", header = TRUE)
sisquoc
#>     name       lon      lat
#> 1    a17 -119.7603 34.75474
#> 2 a20-24 -119.7563 34.75380
#> 3 a25-28 -119.7537 34.75371
#> 4 a18,19 -119.7573 34.75409
#> 5 a35,36 -119.7467 34.75144
#> 6    a31 -119.7478 34.75234
#> 7    a38 -119.7447 34.75230
#> 8    a43 -119.7437 34.75251

# note that ggmap tends to use "lon" instead of "long" for longitude.
ggmap typically asks you for a zoom level, but we can try using ggmap’s make_bbox function:
  
  sbbox <- make_bbox(lon = sisquoc$lon, lat = sisquoc$lat, f = .1)
sbbox
#>       left     bottom      right        top 
#> -119.76198   34.75111 -119.74201   34.75507
Now, when we grab the map ggmap will try to fit it into that bounding box. Let’s try:
  
  # First get the map. By default it gets it from Google.  I want it to be a satellite map
  sq_map <- get_map(location = sbbox, maptype = "satellite", source = "google")
#> Warning: bounding box given to google - spatial extent only approximate.
#> converting bounding box to center/zoom specification. (experimental)
#> Map from URL : http://maps.googleapis.com/maps/api/staticmap?center=34.75309,-119.751995&zoom=16&size=640x640&scale=2&maptype=satellite&language=en-EN&sensor=false
ggmap(sq_map) + geom_point(data = sisquoc, mapping = aes(x = lon, y = lat), color = "red")
#> Warning: Removed 3 rows containing missing values (geom_point).

Nope! That was a fail, but we got a warning about it too. (Actually it is a little better than before because I hacked ggmap a bit…) Let’s try using the zoom level. Zoom levels go from 3 (world scale to 20 (house scale)).

# compute the mean lat and lon
ll_means <- sapply(sisquoc[2:3], mean)
sq_map2 <- get_map(location = ll_means,  maptype = "satellite", source = "google", zoom = 15)
#> Map from URL : http://maps.googleapis.com/maps/api/staticmap?center=34.753117,-119.751324&zoom=15&size=640x640&scale=2&maptype=satellite&language=en-EN&sensor=false
ggmap(sq_map2) + 
  geom_point(data = sisquoc, color = "red", size = 4) +
  geom_text(data = sisquoc, aes(label = paste("  ", as.character(name), sep="")), angle = 60, hjust = 0, color = "yellow")

That is decent. How about if we use the “terrain” type of map:
  
  sq_map3 <- get_map(location = ll_means,  maptype = "terrain", source = "google", zoom = 15)
#> Map from URL : http://maps.googleapis.com/maps/api/staticmap?center=34.753117,-119.751324&zoom=15&size=640x640&scale=2&maptype=terrain&language=en-EN&sensor=false
ggmap(sq_map3) + 
  geom_point(data = sisquoc, color = "red", size = 4) +
  geom_text(data = sisquoc, aes(label = paste("  ", as.character(name), sep="")), angle = 60, hjust = 0, color = "yellow")


That is cool, but I would search for a better color for the lettering…

How about a bike ride?
  I was riding my bike one day with a my phone and downloaded the GPS readings at short intervals.
We can plot the route like this:
  
  bike <- read.csv("data/bike-ride.csv")
head(bike)
#>         lon      lat elevation                 time
#> 1 -122.0646 36.95144      15.8 2011-12-08T19:37:56Z
#> 2 -122.0646 36.95191      15.5 2011-12-08T19:37:59Z
#> 3 -122.0645 36.95201      15.4 2011-12-08T19:38:04Z
#> 4 -122.0645 36.95218      15.5 2011-12-08T19:38:07Z
#> 5 -122.0643 36.95224      15.7 2011-12-08T19:38:10Z
#> 6 -122.0642 36.95233      15.8 2011-12-08T19:38:13Z


bikemap1 <- get_map(location = c(-122.080954, 36.971709), maptype = "terrain", source = "google", zoom = 14)
#> Map from URL : http://maps.googleapis.com/maps/api/staticmap?center=36.971709,-122.080954&zoom=14&size=640x640&scale=2&maptype=terrain&language=en-EN&sensor=false
ggmap(bikemap1) + 
  geom_path(data = bike, aes(color = elevation), size = 3, lineend = "round") + 
  scale_color_gradientn(colours = rainbow(7), breaks = seq(25, 200, by = 25))

See how we have mapped elevation to the color of the path using our rainbow colors again.
Note that getting the right zoom and position for the map is sort of trial and error. You can go to google maps to figure out where the center should be (right click and choose “What’s here?” to get the lat-long of any point. )
The make_bbox function has never really worked for me.

Fish sampling locations
For this, I have whittled down some stuff in the coded wire tag data base to georeferenced marine locations in British Columbia where at least one Chinook salmon was recovered in between 2000 and 2012 inclusive. To see how I did all that you can check out this

Let’s have a look at the data:
  
  bc <- readRDS("data/bc_sites.rds")

# look at some of it:
bc %>% select(state_or_province:sub_location, longitude, latitude)
#> Source: local data frame [1,113 x 9]
#> 
#>    state_or_province water_type sector region area location sub_location
#> 1                  2          M      S     22 016   THOR IS           01
#> 2                  2          M      N     26 012   MITC BY           18
#> 3                  2          M      S     22 015   HARW IS           02
#> 4                  2          M      N     26 006   HOPK PT           01
#> 5                  2          M      S     23 017   TENT IS           06
#> 6                  2          M      S     28 23A   NAHM BY           02
#> 7                  2          M      N     26 006   GIL IS            06
#> 8                  2          M      S     27 024   CLEL IS           06
#> 9                  2          M      S     27 23B   SAND IS           04
#> 10                 2          M      N     26 012   DUVA IS           16
#> ..               ...        ...    ...    ...  ...      ...          ...
#> Variables not shown: longitude (dbl), latitude (dbl)
So, we have 1,113 points to play with.

What do we hope to learn?
  These locations in BC are hierarchically structured. I am basically interested in how close together sites in the same “region” or “area” or “sector” are, and pondering whether it is OK to aggregate fish recoveries at a certain level for the purposes of getting a better overall estimate of the proportion of fish from different hatcheries in these areas.

So, pretty simple stuff. I just want to plot these points on a map, and paint them a different color according to their sector, region, area, etc.
Let’s just enumerate things first, using dplyr:
  
  bc %>% group_by(sector, region, area) %>% tally()
#> Source: local data frame [42 x 4]
#> Groups: sector, region
#> 
#>    sector region area   n
#> 1             48 008    1
#> 2             48 028    1
#> 3             48 311    1
#> 4       N     25 001   33
#> 5       N     25 003   15
#> 6       N     25 004   44
#> 7       N     25 02E    2
#> 8       N     25 02W   34
#> 9       N     26 006   28
#> 10      N     26 007   23
#> 11      N     26 008   10
#> 12      N     26 009   27
#> 13      N     26 010    3
#> 14      N     26 011   11
#> 15      N     26 012   72
#> 16      S     22 013   67
#> 17      S     22 014   58
#> 18      S     22 015   34
#> 19      S     22 016   32
#> 20      S     23 017   53
#> 21      S     23 018   27
#> 22      S     23 028   30
#> 23      S     23 029   14
#> 24      S     23 19A   12
#> 25      S     24 020   44
#> 26      S     24 19B   49
#> 27      S     27 021    4
#> 28      S     27 022    1
#> 29      S     27 023    2
#> 30      S     27 024   38
#> 31      S     27 025   58
#> 32      S     27 026   11
#> 33      S     27 027   23
#> 34      S     27 23B  109
#> 35      S     27 P025   1
#> 36      S     28 23A   23
#> 37      S     61 013   49
#> 38      S     62 014   29
#> 39      S     62 015   17
#> 40      S     62 016   20
#> 41      S     AF 23A    2
#> 42      S     AF SMRV   1
That looks good. It appears like we could probably color code over the whole area down to region, and then down to area within subregions.

Makin’ a map.
Let us try again to use make_bbox() to see if it will work better when used on a large scale.
# compute the bounding box
bc_bbox <- make_bbox(lat = latitude, lon = longitude, data = bc)
bc_bbox
#>       left     bottom      right        top 
#> -133.63297   47.92497 -122.33652   55.80833

# grab the maps from google
bc_big <- get_map(location = bc_bbox, source = "google", maptype = "terrain")
#> Warning: bounding box given to google - spatial extent only approximate.
#> converting bounding box to center/zoom specification. (experimental)
#> Map from URL : http://maps.googleapis.com/maps/api/staticmap?center=51.86665,-127.98475&zoom=6&size=640x640&scale=2&maptype=terrain&language=en-EN&sensor=false

# plot the points and color them by sector
ggmap(bc_big) + 
  geom_point(data = bc, mapping = aes(x = longitude, y = latitude, color = sector))


Cool! That was about as easy as could be. North is in the north, south is in the south, and the three reddish points are clearly aberrant ones at the mouths of rivers.
Coloring it by region
We should be able to color these all by region to some extent (it might get overwhelming), but let us have a go with it.
Notice that region names are unique overall (not just within N or S) so we can just color by region name.
ggmap(bc_big) + 
  geom_point(data = bc, mapping = aes(x = longitude, y = latitude, color = region))


Once again that was dirt easy, though at this scale with all the different regions, it is hard to resolve all the colors.
Zooming in on each region and coloring by area
It is time to really put this thing through its paces. (Keeping in mind that make_bbox() might fail…)
I want to make series of maps. One for each region, in which the the areas in that region are colored differently.
How? Let’s make a function: you pass it the region and it makes the plot.
Keep in mind that there are no factors in this data frame so we don’t have to worry about dropping levels, etc.
region_plot <- function(MyRegion) {
  tmp <- bc %>% filter(region == MyRegion)
  bbox <- make_bbox(lon = longitude, lat = latitude, data = tmp)
  mymap <- get_map(location = bbox, source = "google", maptype = "terrain")
  # now we want to count up how many areas there are
  NumAreas <- tmp %>% summarise(n_distinct(area))
  NumPoints <- nrow(tmp)
  
  the_map <- ggmap(mymap) +
    geom_point(data = tmp, mapping = aes(x = longitude, y = latitude), size = 4, color = "black") +
    geom_point(data = tmp, mapping = aes(x = longitude, y = latitude, color = area), size = 3) +
    ggtitle(
      paste("BC Region: ", MyRegion, " with ", NumPoints, " locations in ", NumAreas, " area(s)", sep = "")
    )
  
  ggsave(paste("bc_region", MyRegion, ".pdf", sep = ""), the_map, width = 9, height = 9)
}
So, with that function we just need to cycle over the regions and make all those plots.

Note that I am saving them to PDFs because it is no fun to make a web page with all of those in there.

dump <- lapply(unique(bc$region), region_plot)



########################################

USA<-map_data("usa")


ggplot() + 
  geom_polygon(data = USA, aes(x=long, y = lat, group = group), fill = NA, color = "black") + 
  coord_fixed(1.3)




#########################################

colr<-c("deeppink","gold",  "purple4","springgreen3","deepskyblue")
plot(net, size=attr(net, "freq")/50, scale.ratio=.5, pie=ind.hap,
     lwd=2, cex=.8, legend=F, fast=F, show.mutation=1, threshold=0, 
     bg=colr)
legend(10, -10, colnames(ind.hap), col=colr, pch=19, cex=1.5)




#also look up diffhaplo
https://stackoverflow.com/questions/49446128/identifying-mutations-between-haplotypes-using-haplonet-pegas-in-r

#pokeweed cpDNA haplotype network
