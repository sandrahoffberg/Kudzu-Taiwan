library(conStruct)

# import allele frequency matrix
crano_freqs <- read.csv("C:/Users/Sandra/OneDrive/kudzu/Taiwan/indiv_TW.fullChina_lessmiss_freqs2.csv", header = TRUE)
crano_freqs <- crano_freqs[, 2:1214]
crano_freqs <- as.matrix(crano_freqs)
#class(crano_freqs)


#import geographic sampling coordinates (longitude, latitude)
crano_coords <- read.csv("C:/Users/Sandra/OneDrive/kudzu/Taiwan/TW.China_3RAD_GPS_lessmiss2.csv", header=T)
crano_coords <- crano_coords[, 1:2]
crano_coords <- as.matrix(crano_coords)
#class(crano_coords)


#create geographic distance matrix
crano_geo <- fields::rdist.earth(crano_coords)
#class(crano_geo)




#### cross-validation ########


# to run a cross-validation analysis
# you have to specify:
#       the numbers of layers you want to compare (K)
#       the allele frequency data (freqs)
#       the geographic distance matrix (geoDist)
#       the sampling coordinates (coords)

setwd("C:/Users/Sandra/OneDrive/kudzu/Taiwan/construct/")
my.xvals <- x.validation(train.prop = 0.9,
                         n.reps = 3,
                         K = 1:10,
                         freqs = crano_freqs,
                         #data.partitions = NULL,
                         geoDist = crano_geo,
                         coords = crano_coords,
                         prefix = "xval_TW.fullChina_unlinked",
                         n.iter = 5000,
                         control = setNames(list(0.99), "adapt_delta"),
                         make.figs = TRUE,
                         save.files = TRUE,
                         parallel = TRUE,
                         n.nodes = 4)



# read in results from text files
#sp.results <- as.matrix(read.table("xval_TW_sp_xval_results.txt",header = TRUE,stringsAsFactors = FALSE))
#nsp.results <- as.matrix(read.table("xval_TW_nsp_xval_results.txt",header = TRUE,stringsAsFactors = FALSE))

# or, format results from the output list
sp.results <- Reduce("cbind",lapply(my.xvals,function(x){unlist(x$sp)}),init=NULL)
nsp.results <- Reduce("cbind",lapply(my.xvals,function(x){unlist(x$nsp)}),init=NULL)




# first, get the 95% confidence intervals for the spatial and nonspatial
#   models over values of K (mean +/- 1.96 the standard error)

sp.CIs <- apply(sp.results,1,function(x){mean(x) + c(-1.96,1.96) * sd(x)/length(x)})
nsp.CIs <- apply(nsp.results,1,function(x){mean(x) + c(-1.96,1.96) * sd(x)/length(x)})

# then, plot cross-validation results for K=1:3 with 8 replicates

par(mfrow=c(1,2))
plot(rowMeans(sp.results),
     pch=19,col="blue",
     ylab="predictive accuracy",xlab="values of K",
     ylim=range(sp.results,nsp.results),
     main="cross-validation results")
points(rowMeans(nsp.results),col="green",pch=19)

# finally, visualize results for the spatial model
#   separately with its confidence interval bars
#
# note that you could do the same with the spatial model, 
#   but the confidence intervals don't really show up 
#   because the differences between predictive accuracies
#   across values of K are so large.

plot(rowMeans(sp.results),
     pch=19,col="blue",
     ylab="predictive accuracy",xlab="values of K",
     ylim=range(sp.CIs),
     main="spatial cross-validation results")
segments(x0 = 1:nrow(sp.results),
         y0 = sp.CIs[1,],
         x1 = 1:nrow(sp.results),
         y1 = sp.CIs[2,],
         col = "blue",lwd=2)














```{r}
#run conStruct analysis
crano <- conStruct(spatial = TRUE,
                   K = 2,
                   freqs = crano_freqs,
                   geoDist = crano_geo,
                   coords = crano_coords,
                   prefix = "test",
                   n.chains = 4,
                   n.iter = 100,
                   make.figs = T,
                   save.files = T,
                   control = setNames(list(0.99), "adapt_delta")
                   # control = setNames(list(15),"max_treedepth")
)
```

###### visualize the results with structure plot  for K=2 #######

load("C:/Users/Sandra/OneDrive/kudzu/Taiwan/construct/TW_K2_10k_data.block.Robj")
load("C:/Users/Sandra/OneDrive/kudzu/Taiwan/construct/TW_K2_10k_construct.results.Robj")
K2_cr <-conStruct.results
K2_db <-data.block


admix.props1 <- conStruct.results$chain_1$MAP$admix.proportions
admix.props2 <- conStruct.results$chain_2$MAP$admix.proportions
admix.props3 <- conStruct.results$chain_3$MAP$admix.proportions
admix.props4 <- conStruct.results$chain_4$MAP$admix.proportions




#match up layers between runs

match.layers.x.runs(admix.props1,admix.props4)



sampnames = c(  "KTW1-11A",  "KTW1-14A",  "KTW1-15A",  "KTW1-1A",  "KTW1-22A",  "KTW1-23A",  "KTW1-3A",  "KTW1-4A",  "KTW1-7A",  "KTW1-9A",  "KTW4-12B",  "KTW4-13B",  "KTW4-15B",  "KTW4-17B",  "KTW4-18B",   "KTW4-19B",   "KTW4-20B",   "KTW4-21B",   "KTW4-22B",   "KTW4-24B",   "KTW6-12A",   "KTW6-15A",   "KTW6-17A",   "KTW6-18A",   "KTW6-1A",   "KTW6-23A",   "KTW6-24A",   "KTW6-2A",   "KTW6-5A",   "KTW6-7A",   "KFU7-10K",   "KFU7-11K",   "KFU7-12K",   "KFU7-16K",   "KFU7-21K",   "KFU7-25K",   "KFU7-5A",   "KFU7-7KA",   "KFU7-8KA",   "KFU7-2KA")

make.structure.plot(admix.proportions = admix.props1, 
                    sample.order = c(1:40), sample.names = samporder,
                    c(4,4,2,2), layer.order = c(1,2))

make.structure.plot(admix.proportions = admix.props2, 
                    sample.order = c(1:40), sample.names = samporder,
                    c(4,4,2,2), layer.order = c(2,1), layer.colors = c("red", "blue"))

make.structure.plot(admix.proportions = admix.props3, 
                    sample.order = c(1:40), sample.names = samporder,
                    c(4,4,2,2), layer.order = c(2,1), layer.colors = c("red", "blue"))

make.structure.plot(admix.proportions = admix.props4, 
                    sample.order = c(1:40), sample.names = samporder,
                    c(4,4,2,2), layer.order = c(1,2))


#visualize the results with admix pie plot 
maps::map(xlim = range(data.block$coords[,1]) + c(-2,1), ylim = range(data.block$coords[,2])+c(-1,1), col="black")

make.admix.pie.plot(admix.proportions = admix.props1, 
                    coords = crano_coords, add = TRUE)


make.admix.pie.plot(admix.proportions = admix.props2, 
                    coords = crano_coords, add = TRUE)


make.admix.pie.plot(admix.proportions = admix.props3, 
                    coords = crano_coords, add = TRUE)

make.admix.pie.plot(admix.proportions = admix.props4, 
                    coords = crano_coords, add = TRUE)




###### visualize the results with structure plot  for K=3 #######

load("C:/Users/Sandra/OneDrive/kudzu/Taiwan/construct/TW_K3_10k_data.block.Robj")
load("C:/Users/Sandra/OneDrive/kudzu/Taiwan/construct/TW_K3_10k_construct.results.Robj")
K3_cr <-conStruct.results
K3_db <-data.block


admix.props.k3.1 <- K3_cr$chain_1$MAP$admix.proportions
admix.props.k3.2 <- K3_cr$chain_2$MAP$admix.proportions
admix.props.k3.3 <- K3_cr$chain_3$MAP$admix.proportions
admix.props.k3.4 <- K3_cr$chain_4$MAP$admix.proportions




#match up layers between runs

match.layers.x.runs(admix.props.k3.1,admix.props.k3.4)

sampnames = c(  "KTW1-11A",  "KTW1-14A",  "KTW1-15A",  "KTW1-1A",  "KTW1-22A",  "KTW1-23A",  "KTW1-3A",  "KTW1-4A",  "KTW1-7A",  "KTW1-9A",  "KTW4-12B",  "KTW4-13B",  "KTW4-15B",  "KTW4-17B",  "KTW4-18B",   "KTW4-19B",   "KTW4-20B",   "KTW4-21B",   "KTW4-22B",   "KTW4-24B",   "KTW6-12A",   "KTW6-15A",   "KTW6-17A",   "KTW6-18A",   "KTW6-1A",   "KTW6-23A",   "KTW6-24A",   "KTW6-2A",   "KTW6-5A",   "KTW6-7A",   "KFU7-10K",   "KFU7-11K",   "KFU7-12K",   "KFU7-16K",   "KFU7-21K",   "KFU7-25K",   "KFU7-5A",   "KFU7-7KA",   "KFU7-8KA",   "KFU7-2KA")

make.structure.plot(admix.proportions = admix.props.k3.1, 
                    sample.order = c(1:40), sample.names = samporder,
                    c(4,4,2,2), layer.order = c(1,2,3), layer.colors = c("blue","red", "gold"))

make.structure.plot(admix.proportions = admix.props.k3.2, 
                    sample.order = c(1:40), sample.names = samporder,
                    c(4,4,2,2), layer.order = c(2,3,1), layer.colors = c("gold","blue", "red"))

make.structure.plot(admix.proportions = admix.props.k3.3, 
                    sample.order = c(1:40), sample.names = samporder,
                    c(4,4,2,2), layer.order = c(2,3,1), layer.colors = c("gold", "blue", "red"))

make.structure.plot(admix.proportions = admix.props.k3.4, 
                    sample.order = c(1:40), sample.names = samporder,
                    c(4,4,2,2), layer.order = c(3,1,2), layer.colors = c("red", "gold","blue"))


#visualize the results with admix pie plot 
maps::map(xlim = range(data.block$coords[,1]) + c(-2,1), ylim = range(data.block$coords[,2])+c(-1,1), col="black")

make.admix.pie.plot(admix.proportions = admix.props.k3.1, 
                    coords = crano_coords, add = TRUE)


make.admix.pie.plot(admix.proportions = admix.props.k3.2, 
                    coords = crano_coords, add = TRUE)


make.admix.pie.plot(admix.proportions = admix.props.k3.3, 
                    coords = crano_coords, add = TRUE)

make.admix.pie.plot(admix.proportions = admix.props.k3.4, 
                    coords = crano_coords, add = TRUE)




###### visualize the results with structure plot  for K=4 #######

load("C:/Users/Sandra/OneDrive/kudzu/Taiwan/construct/TW_K4_10k_data.block.Robj")
load("C:/Users/Sandra/OneDrive/kudzu/Taiwan/construct/TW_K4_10k_construct.results.Robj")
K4_cr <-conStruct.results
K4_db <-data.block


admix.props.k4.1 <- K4_cr$chain_1$MAP$admix.proportions
admix.props.k4.2 <- K4_cr$chain_2$MAP$admix.proportions
admix.props.k4.3 <- K4_cr$chain_3$MAP$admix.proportions
admix.props.k4.4 <- K4_cr$chain_4$MAP$admix.proportions




#match up layers between runs

match.layers.x.runs(admix.props.k4.1,admix.props.k4.4)

sampnames = c(  "KTW1-11A",  "KTW1-14A",  "KTW1-15A",  "KTW1-1A",  "KTW1-22A",  "KTW1-23A",  "KTW1-3A",  "KTW1-4A",  "KTW1-7A",  "KTW1-9A",  "KTW4-12B",  "KTW4-13B",  "KTW4-15B",  "KTW4-17B",  "KTW4-18B",   "KTW4-19B",   "KTW4-20B",   "KTW4-21B",   "KTW4-22B",   "KTW4-24B",   "KTW6-12A",   "KTW6-15A",   "KTW6-17A",   "KTW6-18A",   "KTW6-1A",   "KTW6-23A",   "KTW6-24A",   "KTW6-2A",   "KTW6-5A",   "KTW6-7A",   "KFU7-10K",   "KFU7-11K",   "KFU7-12K",   "KFU7-16K",   "KFU7-21K",   "KFU7-25K",   "KFU7-5A",   "KFU7-7KA",   "KFU7-8KA",   "KFU7-2KA")

make.structure.plot(admix.proportions = admix.props.k4.1, 
                    sample.order = c(1:40), sample.names = samporder,
                    c(4,4,2,2), layer.order = c(1,2,3,4))#, layer.colors = c("blue","red", "gold"))

make.structure.plot(admix.proportions = admix.props.k4.2, 
                    sample.order = c(1:40), sample.names = samporder,
                    c(4,4,2,2), layer.order = c(4,2,1,3), layer.colors = c("red","gold", "forestgreen", "blue"))

make.structure.plot(admix.proportions = admix.props.k4.3, 
                    sample.order = c(1:40), sample.names = samporder,
                    c(4,4,2,2), layer.order = c(1,4,2,3), layer.colors = c("blue", "gold","red", "forestgreen"))

make.structure.plot(admix.proportions = admix.props.k4.4, 
                    sample.order = c(1:40), sample.names = samporder,
                    c(4,4,2,2), layer.order = c(1,4,2,3), layer.colors = c("blue", "forestgreen","red", "gold"))


#visualize the results with admix pie plot 
maps::map(xlim = range(data.block$coords[,1]) + c(-2,1), ylim = range(data.block$coords[,2])+c(-1,1), col="black")

make.admix.pie.plot(admix.proportions = admix.props.k4.1, 
                    coords = crano_coords, add = TRUE)


make.admix.pie.plot(admix.proportions = admix.props.k4.2, 
                    coords = crano_coords, add = TRUE)


make.admix.pie.plot(admix.proportions = admix.props.k4.3, 
                    coords = crano_coords, add = TRUE)

make.admix.pie.plot(admix.proportions = admix.props.k4.4, 
                    coords = crano_coords, add = TRUE)






# compare two runs

setwd("C:/Users/Sandra/OneDrive/kudzu/Taiwan/construct")

compare.two.runs(conStruct.results1=K3_cr,
                 data.block1=K3_db,
                 conStruct.results2=K4_cr,
                 data.block2=K4_db,
                 prefix="TW_10K_K3_vs_K4")

# generates a bunch of pdf figures
```





# Loop through output files generated by conStruct 
#   runs with K=1 through 5 and calculate the 
#   layer contributions for each layer in each run  

layer.contributions <- matrix(NA,nrow=5,ncol=5)

# load the conStruct.results.Robj and data.block.Robj
#   files saved at the end of a conStruct run
load("xval_TW_nsp_rep1K2_construct.results.Robj")
load("xval_TW_nsp_rep1K2_data.block.Robj")

# calculate layer contributions
layer.contributions[,1] <- c(calculate.layer.contribution(conStruct.results[[1]],data.block),rep(0,4))
tmp <- conStruct.results[[1]]$MAP$admix.proportions

for(i in 2:5){
  # load the conStruct.results.Robj and data.block.Robj
  #   files saved at the end of a conStruct run
  load(sprintf("K%s_sp_conStruct.results.Robj",i))
  load(sprintf("K%s_sp_data.block.Robj",i))
  
  # match layers up across runs to keep plotting colors consistent
  #   for the same layers in different runs
  tmp.order <- match.layers.x.runs(tmp,conStruct.results[[1]]$MAP$admix.proportions)  
  
  # calculate layer contributions
  layer.contributions[,i] <- c(calculate.layer.contribution(conStruct.results=conStruct.results[[1]],
                                                            data.block=data.block,
                                                            layer.order=tmp.order),
                               rep(0,5-i))
  tmp <- conStruct.results[[1]]$MAP$admix.proportions[,tmp.order]
}










```{r}
# read in results from text files

sp.results <- as.matrix(
  read.table("lachno_sp_xval_results.txt",
             header = TRUE,
             stringsAsFactors = FALSE)
)
nsp.results <- as.matrix(
  read.table("lachno_nsp_xval_results.txt",
             header = TRUE,
             stringsAsFactors = FALSE)
)

# or, format results from the output list
sp.results <- Reduce("cbind",lapply(my.xvals,function(x){unlist(x$sp)}),init=NULL)
nsp.results <- Reduce("cbind",lapply(my.xvals,function(x){unlist(x$nsp)}),init=NULL)
```





