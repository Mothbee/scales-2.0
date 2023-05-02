setwd("/Users/moth/Downloads/SCALE PROJECT/Analysis/R/scales-2.0")

rm(list=ls())

library(dplyr)
library(phytools)
library(geomorph) # for some reason i need xquartz to be open while running it in R now?
#install.packages("viridis")
library(viridis)
#install.packages("ggcorrplot")
library(ggcorrplot)
#install.packages("igraph")
library(igraph)
library(gprofiler2)
#install.packages("geomorph")



####### Read in data ###########

# scale data
males <- read.csv("males.csv",
                  header = T)
females <- read.csv("females.csv",
                    header = T)

# # habitat/ecology data
# # thx Jill
# hab <- read_excel("/Users/michellesu/Desktop/SCALE PROJECT/Analysis/R/michelle_scales_desktop/Anolis_Species_Data.xlsx")
# hab <- select(hab, Species, ecomorph.DLM, Habitat.type, Open.Shade)




####### Read in and prune tree ######
# Poe et al 2017
tree_full <- read.tree("/Users/moth/Dropbox/MAHLER LAB GENERAL SHARED/MAHLER LAB DATA/ANOLIS PHYLO/trees_cleaned/Poe.et.al.2017.Timetree_cleaned.tre")
tree_spp <- as.data.frame(tree_full$tip.label)


########### male data #########
## male tree ##
data.not.tree <- setdiff(males$species, tree_full$tip.label)
tree.not.data <- setdiff(tree_full$tip.label, males$species)

males_intree <- males[!males$species %in% data.not.tree, ]
howmany_m <- nrow(males_intree)
tree.m.pruned <- drop.tip(tree_full, tree.not.data)

#plot(tree.m.pruned, cex=0.5)

# reordering data for tree

# make row names species
rownames(males_intree) <- males_intree$species # warning: setting row names on tibble is deprecated??

# convert to matrix
males_intree <- as.matrix(males_intree)

# reorder data to match tree
males_intree <- males_intree[tree.m.pruned$tip.label,,drop=F]

#plot(tree.m.pruned, cex=0.5)

########### female data #########
## female tree ##
data.not.tree <- setdiff(females$species, tree_full$tip.label)
tree.not.data <- setdiff(tree_full$tip.label, females$species)

females_intree <- females[!females$species %in% data.not.tree, ]
howmany_f <- (nrow(females_intree))
tree.f.pruned <- drop.tip(tree_full, tree.not.data)

# reordering data for tree

# make row names species
rownames(females_intree) <- females_intree$species # warning: setting row names on tibble is deprecated??

# convert to matrix
females_intree <- as.matrix(females_intree)

# reorder data to match tree
females_intree <- females_intree[tree.f.pruned$tip.label,,drop=F]







####### prepping scale data for analysis #######

# removes rownames but still in order

# remove species column
m_data <- males_intree[,c("svl_mean", "head.count_mean", "midline.count_mean", "dorsal.count_mean", "dorsolateral.count_mean", "chin.count_mean", "ventral.count_mean",
                          "foretoes.count.t_mean", "foretoes.count.b_mean", "hindtoes.count.t_mean", "hindtoes.count.b_mean"), drop=F]

f_data <- females_intree[,c("svl_mean", "head.count_mean", "midline.count_mean", "dorsal.count_mean", "dorsolateral.count_mean", "chin.count_mean", "ventral.count_mean",
                            "foretoes.count.t_mean", "foretoes.count.b_mean", "hindtoes.count.t_mean", "hindtoes.count.b_mean"), drop=F]

class(m_data) <- "numeric" # make all numeric
class(f_data) <- "numeric" # make all numeric



# It's called modtest but i use this for EVERYTHING
m_modtest <- m_data[,c("head.count_mean", "midline.count_mean", "dorsal.count_mean", "dorsolateral.count_mean", "chin.count_mean", "ventral.count_mean",
                       "foretoes.count.t_mean", "foretoes.count.b_mean", "hindtoes.count.t_mean", "hindtoes.count.b_mean"),drop=F]


##### phylo.integration ######
# thanks alex

#cor(m_modtest)
# tomo's advice: look at correlation matrix of the contrasts (PIC)

# You were right; this function does not work without modules as input. However, you can indicate every trait as a module, which has the same effect
integration.allmodules<-phylo.integration(m_modtest, phy= tree.m.pruned, partition.gp = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), iter = 999, seed = "random", print.progress = F)
summary(integration.allmodules)



# making plots
int.toplot <- as.array(integration.allmodules[[2]])

# r-PLS
pls_values <- as.matrix(int.toplot)
pls_values[pls_values==0]<-1

###  heatmap #####
colnames(pls_values) <- c("Snout", "Midline", "Dorsal", "Lateral", "Chin", "Ventral", "Toepad (fore)", "Toe base (fore)", "Toepad (hind)", "Toe base (hind)")
rownames(pls_values) <- c("Snout", "Midline", "Dorsal", "Lateral", "Chin", "Ventral", "Toepad (fore)", "Toe base (fore)", "Toepad (hind)", "Toe base (hind)")

#### ggcorplot
# oh thank god it works
ggcorrplot(
  pls_values,
  type="lower",
  outline.color = "black",
  ggtheme = ggplot2::theme_classic,
  lab_col = "black",
  lab_size = 3.6, # adjust depending on final plot size
  lab=T,
  tl.cex=10,
  tl.col =) +
  scale_fill_viridis()+
  labs(fill="r-PLS")


##### igraph, thanks tomo
# didn't get it to work but maybe one day

# pls_igraph <- pls_values
# pls_igraph[ pls_igraph < .75 ] <- 0
# diag(pls_igraph) <- 0
#
# #labels <- c("Snout", "Midline", "Dorsal", "Lateral", "Chin", "Ventral", "Toepad (fore)", "Toe base (fore)", "Toepad (hind)", "Toe base (hind)")
#
# graph <- graph.adjacency(pls_igraph,add.colnames='label', weighted=TRUE, mode="lower")
# E(graph)$width <- (E(graph)$weight+.5)^8
#
# V(graph)$modules <- c("D", "D", "D", "D","C","C","A","B","A","B")
#
#
# colrs <- c("tomato", "gold", "green", "cyan")
#
# V(graph)$color <- colrs[V(graph)$modules]
#
# plot(graph, layout=layout_nicely(graph))



# #### do all pairwise phylo.integration and use compare.pls ####
# # didn't use
# m_modtest <- m_data[,c("head.count_mean", "midline.count_mean", "dorsal.count_mean", "dorsolateral.count_mean", "chin.count_mean", "ventral.count_mean",
#                        "foretoes.count.t_mean", "foretoes.count.b_mean", "hindtoes.count.t_mean", "hindtoes.count.b_mean"),drop=F]
#
# test1 <- m_data[,c("midline.count_mean", "dorsal.count_mean")]
# test2 <- m_data[,c("hindtoes.count.t_mean", "dorsal.count_mean")]
# test3 <- m_data[,c("ventral.count_mean", "dorsal.count_mean")]
#
#
# int.test1<-phylo.integration(test1, phy= tree.m.pruned, partition.gp = colnames(test1), iter = 999, seed = "random", print.progress = F)
# int.test2<-phylo.integration(test2, phy= tree.m.pruned, partition.gp = colnames(test2), iter = 999, seed = "random", print.progress = F)
# int.test3<-phylo.integration(test3, phy= tree.m.pruned, partition.gp = colnames(test3), iter = 999, seed = "random", print.progress = F)
#
# int.test1
# int.test2
# int.test3
#
# compare.pls(int.test1, int.test2, int.test3) # test 1 diff from test 2


#####  Integration via Eigenvalues from alex ####
#	see Pavlicev et al. (2009) Measuring Morphological Integration Using Eigenvalue Variance
# I think I compute these the way the original author intended, but the resulting integration values are markedly lower than expected.
PCA<-prcomp(m_modtest)
scaledEigen <- PCA$sdev/mean(PCA$sdev)	# Eigenvalues centered around a mean of 1
var(scaledEigen)/(ncol(m_modtest)-1)	# relative Eigenvalue variance, maximum integration at 1






####### phylo.modularity #########
m_modtest <- m_data[,c("head.count_mean", "midline.count_mean", "dorsal.count_mean", "dorsolateral.count_mean", "chin.count_mean", "ventral.count_mean",
                       "foretoes.count.t_mean", "foretoes.count.b_mean", "hindtoes.count.t_mean", "hindtoes.count.b_mean"),drop=F]


####### modularity hypotheses
# use m_modtest and f_modtest


### males ###

### development ###
sink("output/int_dev.txt")

mods <- c("D", "D", "D", "D","C","C","A","B","A","B")

dev.m <- phylo.modularity(A = m_modtest, partition.gp = mods, phy=tree.m.pruned, iter=999)
summary(dev.m)

dev.m[[6]] # specific CRs
dev.nulls <- dev.m[[7]]

sink()


### function ###

## toepads with ventral
sink("output/int_funct1.txt")

mods <- c("A", "A", "A", "A","B","B","B","B","B","B")
funct1.m <- phylo.modularity(A = m_modtest, partition.gp = mods, phy=tree.m.pruned, iter=999, print.progress=F)
summary(funct1.m) # more between than within which actually means more integration??

dev.m[[6]] # specific CRs
funct1.nulls <- funct1.m[[7]]

sink()


## toepads a separate module
sink("output/int_funct2.txt")

mods <- c("A", "A", "A", "A","B","B","C","B","C","B")
funct2.m <- phylo.modularity(A = m_modtest, partition.gp = mods, phy=tree.m.pruned, iter=999, print.progress=F)
summary(funct2.m) # better than below

funct2.m[[6]] # specific CRs
funct2.nulls <- funct2.m[[7]]

sink()


## head, body, limbs
sink("output/int_funct3.txt")

mods <- c("A", "B", "B", "B","A","B","C","C","C","C")
funct3.m <- phylo.modularity(A = m_modtest, partition.gp = mods, phy=tree.m.pruned, iter=999, print.progress=F)
summary(funct3.m)

funct3.m[[6]] # specific CRs
funct3.nulls <- funct3.m[[7]]

sink()



### null histograms
all_nulls <- matrix(c(funct1.nulls, funct2.nulls, funct3.nulls, dev.nulls), ncol=4, dimnames=(list(seq.int(1:1000), c("funct1", "funct2", "funct3", "dev"))))
all_nulls <- as.data.frame(all_nulls)

par(mfrow = c(2,2))
hist(all_nulls$dev ,main="Development",xlab="Permuted Null CR",col="grey",plot = TRUE, xlim=range(1.25,1.5)); segments(x0 = 1.2688, x1 = 1.2688, y0 = 0, y1 = 160, lwd = 2, col = "red"); text(1.24, 180, "1.2688**", cex=0.85, adj=0)
hist(all_nulls$funct1 ,main="Function 1",xlab="Permuted Null CR",col="grey",plot = TRUE, xlim=range(1.04,1.18)); segments(x0 = 1.066, x1 = 1.066, y0 = 0, y1 = 274, lwd = 2, col = "red"); text(1.055, 309, "1.066*", cex=0.85, adj=0)
hist(all_nulls$funct2 ,main="Function 2",xlab="Permuted Null CR",col="grey",plot = TRUE, xlim=range(1.15,1.35)); segments(x0 = 1.1636, x1 = 1.1636, y0 = 0, y1 = 210, lwd = 2, col = "red"); text(1.145, 236, "1.1636**", cex=0.85, adj=0)
hist(all_nulls$funct3 ,main="Function 3",xlab="Permuted Null CR",col="grey",plot = TRUE, xlim=range(1.15,1.35)); segments(x0 = 1.167, x1 = 1.167, y0 = 0, y1 = 140, lwd = 2, col = "red"); text(1.15, 156, "1.167**", cex=0.85, adj=0)



#### compare.CR ####
# more negative effects represents stronger modular signal
# functional 1, 2, 3
sink("output/compare_CR_output.txt")
compare.CR(dev.m, funct1.m, funct2.m, funct3.m, CR.null=F)
sink()



####### compare.multi.evol.rates ########
# the rate ratio for more than 2 subunits is just max rate/min rate out of all subunits

### males ###
# using this again
m_modtest <- m_data[,c("head.count_mean", "midline.count_mean", "dorsal.count_mean", "dorsolateral.count_mean", "chin.count_mean", "ventral.count_mean",
                       "foretoes.count.t_mean", "foretoes.count.b_mean", "hindtoes.count.t_mean", "hindtoes.count.b_mean"),drop=F]


### development ###
sink("output/rates_dev_pairwise_p.txt")

mods <- c("D", "D", "D", "D","C","C","A","B","A","B")
rates.dev.m <- compare.multi.evol.rates(A = m_modtest, gp = mods, phy=tree.m.pruned, Subset=F, iter=999)
summary(rates.dev.m) # p 0.001, dorsal fastest again, then ventral, then toebase, then toepads truly slowest

rates.dev.m[[5]] # pairwise ratios
rates.dev.m[[6]] # pairwise p

sink()

### ecology ###

## toepads with ventral
sink("output/rates_funct1_pairwise_p.txt")

mods <- c("A", "A", "A", "A","B","B","B","B","B","B")
rates.funct1.m <- compare.multi.evol.rates(A = m_modtest, gp = mods, phy=tree.m.pruned, Subset=F, iter=999)
summary(rates.funct1.m) # not sig

rates.funct1.m[[5]] # pairwise ratios
rates.funct1.m[[6]] # pairwise p

sink()


## toepads a separate module (my golden child)
sink("output/rates_funct2.txt")

mods <- c("A", "A", "A", "A","B","B","C","B","C","B")
rates.funct2.m <- compare.multi.evol.rates(A = m_modtest, gp = mods, phy=tree.m.pruned, Subset=F, iter=999)
summary(rates.funct2.m)  # p 0.001, dorsal module fastest, toepads slowest

rates.funct2.m[[5]] # pairwise ratios
rates.funct2.m[[6]] # pairwise p

sink()


# head, body, limbs
sink("output/rates_funct3_pairwise_p.txt")

mods <- c("A", "B", "B", "B","A","B","C","C","C","C")
rates.funct3.m <- compare.multi.evol.rates(A = m_modtest, gp = mods, phy=tree.m.pruned, Subset=F, iter=999)
summary(rates.funct3.m) # not sig

rates.funct3.m[[5]] # pairwise ratios
rates.funct3.m[[6]] # pairwise p

sink()



### everything separate ###

sink("rates_all_pairwise_p.txt")

test_allrates <- compare.multi.evol.rates(A = m_modtest, gp = colnames(m_modtest), phy=tree.m.pruned, Subset=F, iter=999)
summary(test_allrates) # makes sense

test_allrates[[5]] # pairwise ratios
test_allrates[[6]] # pairwise p

sink()


#### ggcorplot but it's for rate ratios
# maybe not the best way to represent this
all_rateratio <- as.matrix(test_allrates[[5]])
all_rateratio[all_rateratio==0] <- 1
all_rateratio_ps <- as.matrix(test_allrates[[6]])

colnames(all_rateratio) <- c("Snout", "Midline", "Dorsal", "Lateral", "Chin", "Ventral", "Toepad (fore)", "Toe base (fore)", "Toepad (hind)", "Toe base (hind)")
rownames(all_rateratio) <- c("Snout", "Midline", "Dorsal", "Lateral", "Chin", "Ventral", "Toepad (fore)", "Toe base (fore)", "Toepad (hind)", "Toe base (hind)")

ggcorrplot(
  all_rateratio,
  type="lower",
  outline.color = "white",
  #ggtheme = ggplot2::theme_bw,
  p.mat = all_rateratio_ps,
  sig.level = 0.05,
  insig = "blank",
  #pch=4,
  pch.cex =25,
  lab_col = "black",
  lab_size = 7, # adjust depending on final plot size
  lab=T) +
  scale_fill_viridis()+
  labs(fill="Rate ratio")





##### getting colours for heatmap ####
all_rates <- test_allrates[[4]]

#map viridis to numbers
vir_map <- mapViridis(all_rates, domain_min = max(all_rates), domain_max = min(all_rates), n=256)
vir_map
mapped = c("#25848EFF", "#FDE725FF", "#D8E219FF", "#31678EFF", "#470E61FF", "#20A486FF", "#404688FF", "#440154FF", "#A8DB34FF", "#34608DFF")

# just to plot the gradient
imgoingcrazy <- data.frame(x=as.factor(1:10), y= all_rates)
colnames(imgoingcrazy) <- c("region", "Rate")

t<-ggplot(imgoingcrazy, aes(region, Rate, fill=Rate))
t + geom_bar(stat = "identity") +scale_fill_viridis()



### same for dev
dev.rates <- rates.dev.m[[4]]

#map viridis to numbers
vir_map <- mapViridis(dev.rates, domain_min = max(all_rates), domain_max = min(all_rates), n=256)

# just to plot the gradient
imgoingcrazy <- data.frame(x=as.factor(1:4), y= dev.rates)
colnames(imgoingcrazy) <- c("region", "Rate")

t<-ggplot(imgoingcrazy, aes(region, Rate, fill=Rate))
t + geom_bar(stat = "identity") +scale_fill_viridis()







###### thank you alex ######
### Alex Tinius, 6.8.2021; RateFrameOld converts the output of "compare.multi.evol.rates" into a data frame
RateFrameOld <- function(DisInputFull, ROUND=3) {	# ROUND indicates number of decimals retained
  DisInput1<-DisInputFull$sigma.d.gp.ratio
  DisInput2<-DisInputFull$pairwise.pvalue
  DisLabel<-DisInputFull$groups

  DisTable<-NULL
  maxLabel<-length(DisLabel)-1
  for(i in 1:maxLabel){
    if(i==1){DisFirst<-1}
    DisCells<-c(DisFirst:{DisFirst+maxLabel-i})
    DisFirst<-DisFirst+length(DisCells)
    DisCol<-c(rep("",{maxLabel-length(DisCells)}),round(DisInput1[DisCells],ROUND))
    DisTable<-cbind(DisTable,DisCol)
  }
  DisTable<-data.frame(cbind(rbind(rep("",maxLabel),DisTable),rep("",maxLabel+1)))
  colnames(DisTable)<-rownames(DisTable)<-DisLabel

  for(i in 1:maxLabel){
    if(i==1){DisFirst<-1}
    DisCells<-c(DisFirst:{DisFirst+maxLabel-i})	# number of cells with data
    DisFirst<-DisFirst+length(DisCells)	# first cell with data
    DisCol<-c(round(DisInput2[DisCells],ROUND))
    DisTable[i,c((i+1):(i+length(DisCells)))]<-DisCol
  }
  return(DisTable)
}	# end "RateFrameOld" function


### Alex Tinius, 23.8.2021; RateFrame converts the output of "compare.multi.evol.rates" into a data frame, and reorders the matrix
RateFrame <- function(DisInputFull, ROUND=3) {	# ROUND indicates number of decimals retained
  DisInput1<-DisInputFull$sigma.d.gp.ratio	# net evolutionary rate for each group
  DisInput2<-DisInputFull$pairwise.pvalue		# significance level of the observed ratio
  DisLabel<-DisInputFull$groups			# group names
  DisLabelNew<-c(read.table(file="input/_DistLabelNew.txt", sep="\t", dec=".", header=F, row.names=1, skip=1)[,1])	# predefined vector to order the matrix by

  DisTable1<-DisTable2<-DisTable<-NULL
  maxLabel<-length(DisLabel)-1	# since the diagonal elements will be eliminated, the resulting matrix is one column/row shorter
  for(i in 1:maxLabel){
    if(i==1){DisFirst<-1}
    DisCells<-c(DisFirst:{DisFirst+maxLabel-i})
    DisFirst<-DisFirst+length(DisCells)	# first cell containing data in this row
    DisCol<-c(rep("",{maxLabel-length(DisCells)}),round(DisInput1[DisCells],ROUND))	# data for row[i] of the new matrix
    DisTable1<-cbind(DisTable1,DisCol)
  }
  DisTable1<-data.frame(cbind(rbind(rep("",maxLabel),DisTable1),rep("",maxLabel+1)))	# the lower left half of the matrix remains empty
  colnames(DisTable1)<-rownames(DisTable1)<-DisLabel	# rename columns to fit the labels

  for(i in 1:maxLabel){
    if(i==1){DisFirst<-1}
    DisCells<-c(DisFirst:{DisFirst+maxLabel-i})	# number of cells with data
    DisFirst<-DisFirst+length(DisCells)	# first cell with data
    DisCol<-c(round(DisInput1[DisCells],ROUND))
    DisTable1[i,c((i+1):(i+length(DisCells)))]<-DisCol
  }
  DisTable1new<-DisTable1[DisLabelNew,]
  DisTable1new<-DisTable1new[,DisLabelNew]

  # repeat process for second half of the matrix
  for(i in 1:maxLabel){
    if(i==1){DisFirst<-1}
    DisCells<-c(DisFirst:{DisFirst+maxLabel-i})
    DisFirst<-DisFirst+length(DisCells)
    DisCol<-c(rep("",{maxLabel-length(DisCells)}),round(DisInput2[DisCells],ROUND))
    DisTable2<-cbind(DisTable2,DisCol)
  }
  DisTable2<-data.frame(cbind(rbind(rep("",maxLabel),DisTable2),rep("",maxLabel+1)))
  colnames(DisTable2)<-rownames(DisTable2)<-DisLabel

  for(i in 1:maxLabel){
    if(i==1){DisFirst<-1}
    DisCells<-c(DisFirst:{DisFirst+maxLabel-i})	# number of cells with data
    DisFirst<-DisFirst+length(DisCells)	# first cell with data
    DisCol<-c(round(DisInput2[DisCells],ROUND))
    DisTable2[i,c((i+1):(i+length(DisCells)))]<-DisCol
  }
  DisTable2new<-DisTable2[DisLabelNew,]	# order both tables to fit one label vector
  DisTable2new<-DisTable2new[,DisLabelNew]

  # unite the two tables
  DisTable<-DisTable1new
  maxLabel<-nrow(DisTable)
  for(i in 1:maxLabel){
    DisTable[i,i:maxLabel]<-DisTable2new[i,i:maxLabel]	# replace empty elements of Table 1 with contents of Table 2
  }
  return(DisTable)
}	# end "RateFrame" function


### Alex Tinius, 15.11.2022, evolRateMV' calculates evolutionary rate from a Procrustes fit
evolRateMV <- function(TraitMatrix, phyTree, eMod="BM") {	# eMod can be "OU" "BM" or "EB"
  library("mvMORPH")

  eMod<-"BM" # evolutionary model, can be "OU" "BM" or "EB"
  #phyTree<- your tree in phylo format
  #TraitMatrix<- traits in your examined module

  PCA<-prcomp(TraitMatrix)
  procPC<-as.matrix(PCA$x)
  procPC<-procPC[phyTree$tip.label,]
  Dframe<-data.frame(procPC)
  PzModel<-mvgls(procPC~1,Dframe,phyTree,model=eMod,penalty="RidgeArch")	# calculate mvgls model
  return(sum(diag(PzModel$sigma$S)))	# report the sum of the covariance diagonals
}	# end of "evolRateMV" function


#################################
# rates on bottom, p-values on top

##### just dev
dev_rfo <- as.matrix(RateFrameOld(rates.dev.m))
class(dev_rfo) <- "numeric"
colnames(dev_rfo) <- c("Toepads", "Limbs", "Ventral", "Dorsal")
rownames(dev_rfo) <- c("Toepads", "Limbs", "Ventral", "Dorsal")

#### ggcorplot but it's for rate ratios
ggcorrplot(
  dev_rfo,
  type="full",
  outline.color = "white",
  #ggtheme = ggplot2::theme_bw,
  #p.mat = all_rateratio_ps,
  #sig.level = 0.05,
  #insig = "blank",
  #pch=4,
  pch.cex =25,
  lab_col = "black",
  lab_size = 7, # adjust depending on final plot size
  lab=T) +
  scale_fill_viridis()+
  labs(fill="Rate ratio")


##### all traits
allrates_rfo <- as.matrix(RateFrameOld(test_allrates))
class(allrates_rfo) <- "numeric"
colnames(allrates_rfo) <- c("Snout", "Midline", "Dorsal", "Lateral", "Chin", "Ventral", "Toepad (fore)", "Toe base (fore)", "Toepad (hind)", "Toe base (hind)")
rownames(allrates_rfo) <- c("Snout", "Midline", "Dorsal", "Lateral", "Chin", "Ventral", "Toepad (fore)", "Toe base (fore)", "Toepad (hind)", "Toe base (hind)")

#### ggcorplot but it's for rate ratios
ggcorrplot(
  allrates_rfo,
  type="full",
  outline.color = "white",
  #ggtheme = ggplot2::theme_bw,
  #p.mat = all_rateratio_ps,
  #sig.level = 0.05,
  #insig = "blank",
  #pch=4,
  pch.cex =25,
  lab_col = "black",
  lab_size = 7, # adjust depending on final plot size
  lab=T) +
  scale_fill_viridis()+
  labs(fill="Rate ratio")








###### habitat stuff, make.simmap and evolvcv.lite (obsolete?) ######
# compare the rate of dorsal and ventral scale evolution in open vs shade
# most are 0 :( should "both" be grouped with "shade" or "open"?

# merge scale and habitat data
males_intree_hab <- merge(as.data.frame(males_intree), hab, by.x = "species",
                          by.y = "Species", all.x = TRUE, all.y = FALSE, sort=F)
rownames(males_intree_hab) <- males_intree_hab$species # warning: setting row names on tibble is deprecated??

# filter to species w open/shade data
males_intree_hab_filter <- males_intree_hab %>%
  filter(Open.Shade=="shade"|Open.Shade=="open")

### make tree for open/shade dataset
data.not.tree <- setdiff(males_intree_hab_filter$species, tree_full$tip.label)
tree.not.data <- setdiff(tree_full$tip.label, males_intree_hab_filter$species)

males.data.o_s <- males_intree_hab_filter[!males_intree_hab_filter$species %in% data.not.tree, ]
tree.m.pruned.o_s <- drop.tip(tree_full, tree.not.data)

# reordering data for tree

# make row names species
rownames(males.data.o_s) <- males.data.o_s$species # warning: setting row names on tibble is deprecated??

# convert to matrix
males.data.o_s <- as.matrix(males.data.o_s)

# reorder data to match tree
males.data.o_s <- males.data.o_s[tree.m.pruned.o_s$tip.label,,drop=F]


# make vector of tip states
open_shade <- as.data.frame(males.data.o_s)$Open.Shade
names(open_shade) <- as.data.frame(males.data.o_s)$species

# do the test
test <- make.simmap(tree.m.pruned.o_s, x=open_shade)

plot(test)   # this does not look good



