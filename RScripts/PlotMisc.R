## This script is used to plot miscellaneous figures for the simulation. It now includes plots of mean trait trajectory, trait/genetic variance trajectory, allele frequency trajectory, and distribution of D-value on the chromosome. More plot will be added in the future. 
## This script is meant to run with interactive R sessions.

rm(list=ls(all=TRUE))
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  # Multiple plot function
  #
  # ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
  # - cols:   Number of columns in layout
  # - layout: A matrix specifying the layout. If present, 'cols' is ignored.
  #
  # If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
  # then plot 1 will go in the upper left, 2 will go in the upper right, and
  # 3 will go all the way across the bottom.
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
library(ggplot2)
library(reshape)

## Plot trait charateristics across all replicates in the standard model

pooled_trait_mean <- matrix(NA, ncol=10, nrow=200)
pooled_trait_variance <- matrix(NA, ncol=10, nrow=200)
for (k in 1:2) {
  Direction <- c("Plus", "Minus")[k]
  for (j in 1:100){
    setwd(paste0("/fs/cbsubscb10/storage/rl683/TemporalScan/Simulations/NQTL100/SimRep", j, "/ExpRep", Direction, "1"))
    full_trait <- matrix(NA,ncol = 10, nrow = 1000)
    selected_trait_mean <- vector("numeric", length = 10)
    for (i in 1:10){
      assign(paste0("Gen", i,"_Trait"), read.table(paste0("Gen", i, "_Trait.txt"), fill=T))
      full_trait[,i] <- get(paste0("Gen", i,"_Trait"))[,1]
    }  
    
    ## get summary stats from all individuals per generation
    trait_mean <- apply(full_trait, 2, mean)
    trait_variance <- apply(full_trait, 2, var)
    
    ## enter these into matrices
    pooled_trait_mean[j+(k-1)*100,] <- trait_mean
    pooled_trait_variance[j+(k-1)*100,] <- trait_variance
  }
}

## Average trait
pooled_trait_mean<-as.data.frame(pooled_trait_mean)
colnames(pooled_trait_mean)<-1:10
pooled_trait_mean<-cbind(pooled_trait_mean, ReplicateNumber=c(1:200))
pooled_trait_mean<-cbind(pooled_trait_mean, Direction=c(rep("Plus",100), rep("Minus",100)))
pooled_trait_mean <- melt(pooled_trait_mean, id.vars=c("ReplicateNumber", "Direction"), value.name="AverageTraitValue", variable.name="Timepoint")
colnames(pooled_trait_mean)[3:4]<-c("Timepoint", "AverageTraitValue")
p1 <-ggplot(data=pooled_trait_mean) +
  geom_line(aes(x=Timepoint, y=AverageTraitValue, group=ReplicateNumber, color=Direction), alpha=0.6) +
  theme(legend.position="none") +
  xlab("Generation Number") +
  ylab("Average Phenotype Value") +
  theme(axis.title = element_text(size = 30)) +
  theme(axis.text=element_text(size=30)) +
  theme(text = element_text(size=30)) + 
  scale_color_manual(values=c("darkorange2", "purple"))

  


## Trait Variance
pooled_trait_variance<-as.data.frame(pooled_trait_variance)
colnames(pooled_trait_variance)<-1:10
pooled_trait_variance<-cbind(pooled_trait_variance, ReplicateNumber=c(1:200))
pooled_trait_variance<-cbind(pooled_trait_variance, Direction=c(rep("Plus",100), rep("Minus",100)))
pooled_trait_variance <- melt(pooled_trait_variance, id.vars=c("ReplicateNumber", "Direction"), value.name="TraitVariance", variable.name="Timepoint")
colnames(pooled_trait_variance)[3:4]<-c("Timepoint", "TraitVariance")
p2 <-ggplot(data=pooled_trait_variance) +
  geom_line(aes(x=Timepoint, y=TraitVariance, group=ReplicateNumber, color=Direction), alpha=0.6) +
  theme(legend.position="none") +
  xlab("Generation Number") +
  ylab("Phenotypic/Genetic Variance") +
  theme(axis.title = element_text(size = 30)) +
  theme(axis.text=element_text(size=30)) +
  theme(text = element_text(size=30)) +
  scale_color_manual(values=c("darkorange2", "purple"))



# Plotting
setwd("/fs/cbsubscb10/storage/rl683/TemporalScan/Figures/Misc/")
png("TraitMean_NQTL100.png", width = 800, height = 750, units = "px", pointsize = 20)
print(p1)
dev.off()

setwd("/fs/cbsubscb10/storage/rl683/TemporalScan/Figures/Misc/")
png("TraitVariance_NQTL100.png", width = 800, height = 750, units = "px", pointsize = 20)
print(p2)
dev.off()

#########################################################################################
## Plot allele frequency trajectory in one single experiment
library(ggplot2)
library(reshape)
for (i in 1){
  setwd(paste0("/fs/cbsubscb10/storage/rl683/TemporalScan/Simulations/NQTL100/SimRep", i, "/ExpRepPlus1"))
  df1 <- read.table("SampleFreqAllTimepoints.txt", sep = ",", header = T)
  setwd(paste0("/fs/cbsubscb10/storage/rl683/TemporalScan/Simulations/NQTL100/SimRep", i, "/ExpRepMinus1"))
  df2 <- read.table("SampleFreqAllTimepoints.txt", sep = ",", header = T)
  
  colnames(df1)[7:16]<-0:9
  df1<- df1[,c(1, 4, 7:16)]
  df1 <- melt(df1, id.vars=c("ID_perm", "es"), value.name="Frequency", variable.name="Timepoint")
  
  colnames(df2)[7:16]<--9:0
  df2[,7:16]<-df2[,16:7]
  df2<- df2[,c(1, 4, 7:16)]
  df2 <- melt(df2, id.vars=c("ID_perm", "es"), value.name="Frequency", variable.name="Timepoint")
  df <- rbind(df2, df1)
  df[,4]<-df[,4]/100
  colnames(df) <- c("PermanentID", "EffectSize", "Timepoint", "Frequency")
  
  p1 <- ggplot() +
    geom_line(data=df[which(df[,2]==0),], aes(x=Timepoint, y=Frequency, group=PermanentID), colour = "Black", size=0.3,alpha=0.012) +
    geom_line(data=df[which(df[,2]>0),], aes(x=Timepoint, y=Frequency, group=PermanentID), colour = "Red", size=1, alpha=1) +
    geom_line(data=df[which(df[,2]<0),], aes(x=Timepoint, y=Frequency, group=PermanentID), colour = "Blue", size=1, alpha=1) +
    geom_vline(xintercept = 10, size=2) +
    geom_segment(aes(x=11, xend=15, y=0.45, yend=0.45), size=2, arrow = arrow(length = unit(0.5, "cm"))) +
    geom_segment(aes(x=9, xend=5, y=0.45, yend=0.45), size=2, arrow = arrow(length = unit(0.5, "cm"))) +
    geom_label(aes(label="Selection for \n LARGER trait values", x = 14, y = 0.55), size=8, colour='black', alpha=0.7, fontface = "bold") +
    geom_label(aes(label="Selection for \n SMALLER trait values", x = 6, y = 0.55), size=8, colour='black', alpha=0.7, fontface = "bold") +
    scale_x_discrete(breaks=c(-9:9), labels=c(seq(10,1,-1), seq(2,10,1)), name="Generation Number") +
    scale_y_continuous(breaks=seq(0,1,0.25), labels=seq(0,1,0.25), name="Allele Frequency") +
    theme(axis.text = element_text(size=30), axis.title = element_text(size=30)) +
    theme(legend.position="none") +
    theme(text = element_text(size=30))
    
  
  setwd("/fs/cbsubscb10/storage/rl683/TemporalScan/Figures/Misc/")
  png(paste0("FrequencyTrajectory_NQT100_", i, ".png"), width = 800, height = 750, units = "px", pointsize = 20)
  print(p1)
  dev.off()
}

#########################################################################################
## Plot D distribution in one single experiment
library(ggplot2)
for (i in 1){
  setwd(paste0("/fs/cbsubscb10/storage/rl683/TemporalScan/Simulations/NQTL100/SimRep", i, "/ExpRepPlus1"))
  df1 <- read.table("SampleFreqAllTimepoints.txt", sep = ",", header = T)
  setwd(paste0("/fs/cbsubscb10/storage/rl683/TemporalScan/Simulations/NQTL100/SimRep", i, "/ExpRepMinus1"))
  df2 <- read.table("SampleFreqAllTimepoints.txt", sep = ",", header = T)
  
  df1 <- df1[,c(1, 3, 4, 16)]
  colnames(df1) <- c("PermanentID", "Postion", "EffectSize", "FrequencyPlus")
  df2 <- df2[,c(1, 3, 4, 16)]
  colnames(df2) <- c("PermanentID", "Position", "EffectSize", "FrequencyMinus")
  df <- merge(df1, df2, by="PermanentID", all=T)
  Position<-apply(df[,c(2,5)], 1, function(x){mean(x, na.rm = T)})
  EffectSize<-apply(df[,c(3,6)], 1, function(x){mean(x, na.rm = T)})
  df <- cbind(df[,c(1,4,7)], Position, EffectSize)
  df[is.na(df)]<-0
  df <- cbind(df, D=abs(df$FrequencyPlus-df$FrequencyMinus)) 
  
  p1 <- ggplot() +
    geom_point(data=df[which(df[,5]==0),], aes(x=Position/1000000, y=D/100), colour = "Black", size=2, alpha=0.1) +
    geom_point(data=df[which(df[,5]>0),], aes(x=Position/1000000, y=D/100), colour = "Red", size=5) +
    geom_point(data=df[which(df[,5]<0),], aes(x=Position/1000000, y=D/100), colour = "Blue", size=5) +
    theme(legend.position="none") +
    xlab("Position on Chromosome (Mbp)") +
    ylab("D-value") +
    theme(axis.title = element_text(size = 30)) +
    theme(axis.text=element_text(size=30)) +
    theme(text = element_text(size=30))
    
  p2 <- ggplot() +
    geom_point(data=df[which(df[,5]==0),], aes(x=FrequencyPlus, y=FrequencyMinus), colour = "Black", alpha=0.1) +
    geom_point(data=df[which(df[,5]>0),], aes(x=FrequencyPlus, y=FrequencyMinus), colour = "Red", size=3) +
    geom_point(data=df[which(df[,5]<0),], aes(x=FrequencyPlus, y=FrequencyMinus), colour = "Blue", size=3) +
    theme(legend.position="none")
  
  setwd("/fs/cbsubscb10/storage/rl683/TemporalScan/Figures/Misc/")
  png(paste0("D_NQTL100_", i, ".png"), width = 800, height = 750, units = "px", pointsize = 20)
  print(p1)
  dev.off()
  
  # setwd("/fs/cbsubscb10/storage/rl683/TemporalScan/Figures/Misc/")
  # png(paste0("EndFrequency_NQTL100_", i, ".png"), width = 800, height = 750, units = "px", pointsize = 20)
  # print(p2)
  # dev.off()
}


####################################################################
## Plot LD decay

position <-seq(0,3e7, length.out =3001)
p<-4*1000*1e-8*position
LD<-(10+p)/(p^2+13*p+22)
plot(LD~position)
library(ggplot2)
df<- data.frame(position, LD)
ggplot(data=df, aes(x=position, y=LD))
ggplot(data=df, aes(x=position, y=LD)) +
  geom_line() +
  scale_x_continuous(limits = c(0,31000000), breaks = seq(0,3e7,length.out =16)) 
ggplot(data=df, aes(x=position, y=LD)) +
  geom_line() +
  scale_x_continuous(limits = c(0,2000000), breaks = seq(0,2000000,length.out =11)) 


junk<-matrix(c(1,2,3,2,3,4,3,4,5,4,5,6,5,6,7,6,7,8), ncol=3, byrow=T)
junk<-as.data.frame(t(junk))
colnames(junk)<-1:6
junk<-cbind(junk, timepoints=1:3)
junk<-melt(junk, id="timepoints")
ggplot(data=junk) +
  geom_line(aes(x=timepoints, y=value, color=variable))
