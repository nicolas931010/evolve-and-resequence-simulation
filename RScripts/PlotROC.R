## This script is used to create ROC curve comparisons among different trait architecture or experimental design scenarios.
## First define the type of comparisons that you want to make. Then, under the following if statement, define the appropriate variables for type of your comparison. Also define the variables outside of the if statement. Startin from line 75, also find the appropriate if statement and follow instructions there. 
## This step is usually very quick, so an interactive R session will be more suitable.
## You might find it easier to write your own R code for your specific purpose to plot the ROC curve. 


rm(list=ls(all=TRUE))

Type <- "QTL_architecture" # can be "QTL_architecture", "WindowSize", "ComputationalMethod", "Replication", "Direction", "Experimental_design"

if (Type == "Direction"){
  NsExpRep <- c(1)
  Directions <- c("Plus", "Minus")
} else if (Type == "Replication"){
  NsExpRep <- c(1,2,5,10)
  Directions <- c("Plus")
} else {
  NExpRep <- 1
  Directions <- c("Plus"
                  , "Minus"
  )
}
WindowSize <- 0
ComputationalMethod <- "D" # comment out if doing method comparison
filename <- "NQTL10_FewerSNPs" # When Type=="Replication", "Direction", "ComputationalMethod"
title <- paste0(Type, "_", filename) # When Type=="Replication", "Direction", "ComputationalMethod"
OutPath="/fs/cbsubscb10/storage/rl683/TemporalScan/" 
SampleSize=50 # make sure that the simulation had SampleSize LARGER or EQUAL to this
NSimRep=100 # make sure that the simulation had NSimRep LARGER or EQUAL to this

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

if (Type %in% c("QTL_architecture", "experimental_design")){ # define filenames, plot title, and plot labels in the begining of the following loop 
  for (i in c(2)){
    titles <- c("Number_of_SNPs"
                , "Number_of_QTLs"
                , "Starting_frequency_10_QTLs"
                , "Starting_frequency_100_QTLs"
                , "Effect_size_distribution_10_QTLs"
                , "Effect_size_distribution_100_QTLs"
                , "Dominance_10_QTLs"
                , "Dominance_100_QTLs"
                , "Position_10_QTLs"
                , "Position_100_QTLs"
                , "Epistasis_scenarios"
    )
    
    filenames <-   list(c("NQTL10", "NQTL10_FewerSNPs")
                        , c("NQTL10", "NQTL2", "NQTL20", "NQTL100", "NQTL200")
                        , c("NQTL10", "NQTL10_FreqLow5", "NQTL10_FreqHigh5")
                        , c("NQTL100", "NQTL100_FreqLow5", "NQTL100_FreqHigh5")
                        , c("NQTL10", "NQTL10_ESDistE")
                        , c("NQTL100", "NQTL100_ESDistE")
                        , c("NQTL10",  "NQTL10_D0", "NQTL10_D1")
                        , c("NQTL100", "NQTL100_D0", "NQTL100_D1")
                        , c("NQTL10", "NQTL10_Clustered")
                        , c("NQTL100", "NQTL100_Clustered")
                        , c("NQTL10", "NQTL10_EpiSce1", "NQTL10_EpiSce2", "NQTL10_EpiSce3", "NQTL10_EpiSce4", "NQTL10_EpiSce5", "NQTL10_EpiSce6")
    )
    
    labels <- list(c("~14000 SNPs (Baseline)", "~1400 SNPs")
                   , c("10 QTLs (Baseline)", "2 QTLs", "20 QTLs", "100 QTLs (Baseline)", "200 QTLs")
                   , c("Random (Baseline)", "Lower than or equal to 5%", "Higher than or equal to 5%")
                   , c("Random (Baseline)", "Lower than or equal to 5%", "Higher than or equal to 5%")
                   , c("Fixed (Baseline)", "Exponential")
                   , c("Fixed (Baseline)", "Exponential")
                   , c("Codominance (Baseline)", "Wild Type in complete dominance", "Mutant in complete dominance")
                   , c("Codominance (Baseline)", "Wild Type in complete dominance", "Mutant in complete dominance")
                   , c("Random (Basline)", "Clustered")
                   , c("Random (Basline)", "Clustered")
                   , c("No epistatisis (Baseline)", "Epistasis scenario 1", "Epistasis scenario 2","Epistasis scenario 3","Epistasis scenario 4","Epistasis scenario 5","Epistasis scenario 6")
    )
    
    title <- titles[i]
    filenames <- filenames[[i]]
    labels <- labels[[i]]
    
    
    n_files <- length(filenames)
    for (Direction in Directions) {
      for (filename in filenames){
        for (k in 1:NSimRep){
          if (Direction == "Plus") {
            DirectionName <- "SingleDirection"
          } else if (Direction == "Minus"){
            DirectionName <- "OppositeDirections"
          }
          setwd(paste0(OutPath, "Simulations/", filename, "/SimRep", k, "/ExpRep", Direction, NExpRep))
          ROC <- read.table("ROC_TwoTimepoints.txt", sep = ",",  header = T, stringsAsFactors = F)
          colnames(ROC)<- c("Proportion_of_genetic_variance_in_gen_1", 
                            "Proportion_of_QTL_detected_weighted", 
                            "False_positive_rate",
                            "Window_size",
                            "Computational_method",
                            "Numbers_of_timepoints_sampled")
          ROC <- ROC[which(ROC$Window_size==WindowSize & ROC$Computational_method==ComputationalMethod),]
          if (k == 1){
            ROC_final <- ROC[,1:3]
          } else {
            ROC_final <- ROC[,1:3]+ROC_final
          }
        }
        ROC_final <- ROC_final/NSimRep
        ROC_final <- cbind(ROC_final, filename)
        color <- which(filenames == filename)
        ROC_final <- cbind(ROC_final, color) 
        if (!exists("ROC_QTL_archetecture_pooled")){
          ROC_QTL_archetecture_pooled <- ROC_final
        } else {
          ROC_QTL_archetecture_pooled <- rbind(ROC_QTL_archetecture_pooled, ROC_final)
        } 
      }
      library(ggplot2)
      
      
      p1 <- ggplot(ROC_QTL_archetecture_pooled, aes(False_positive_rate, Proportion_of_genetic_variance_in_gen_1, color = factor(color, labels = labels))) + 
        #geom_point(size = 3, alpha = 0.8, shape=1) + 
        geom_line(size = 2, alpha = 0.8, linetype="dotted") + 
        #geom_point(data=ROC_QTL_archetecture_pooled[which(ROC_QTL_archetecture_pooled[,1]>0 | ROC_QTL_archetecture_pooled[,3]>0),], aes(False_positive_rate, Proportion_of_genetic_variance_in_gen_1, color = factor(color, labels = labels)), size = 3, alpha = 0.8) + 
        geom_line(data=ROC_QTL_archetecture_pooled[which(ROC_QTL_archetecture_pooled[,1]>0 | ROC_QTL_archetecture_pooled[,3]>0),], aes(False_positive_rate, Proportion_of_genetic_variance_in_gen_1, color = factor(color, labels = labels)), size = 2, alpha = 1) + 
        scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) + 
        xlab("False positive rate") +
        ylab("Proportion of genetic variance in the first generation") +
        theme(legend.position="none") +
        theme(axis.title = element_text(size = 30)) +
        theme(axis.text=element_text(size=30)) +
        theme(text = element_text(size=30)) +
        ggtitle(NULL)
      
      p2 <- ggplot(ROC_QTL_archetecture_pooled, aes(False_positive_rate, Proportion_of_genetic_variance_in_gen_1, color = factor(color, labels = labels))) + 
        #geom_point(size = 3, alpha = 0.8, shape=1) + 
        geom_line(size = 2, alpha = 0.8, linetype="dotted") + 
        #geom_point(data=ROC_QTL_archetecture_pooled[which(ROC_QTL_archetecture_pooled[,1]>0 | ROC_QTL_archetecture_pooled[,3]>0),], aes(False_positive_rate, Proportion_of_genetic_variance_in_gen_1, color = factor(color, labels = labels)), size = 3, alpha = 0.8) + 
        geom_line(data=ROC_QTL_archetecture_pooled[which(ROC_QTL_archetecture_pooled[,1]>0 | ROC_QTL_archetecture_pooled[,3]>0),], aes(False_positive_rate, Proportion_of_genetic_variance_in_gen_1, color = factor(color, labels = labels)), size = 2, alpha = 1) + 
        scale_x_continuous(limits = c(0,0.105), breaks = seq(0,0.1,0.02)) + 
        #scale_y_continuous(limits = c(0,1)) + 
        xlab("False positive rate") +
        ylab("Proportion of genetic variance in the first generation") +
        theme(legend.position="none") +
        theme(axis.title = element_text(size = 30)) +
        theme(axis.text=element_text(size=30)) +
        theme(text = element_text(size=30)) +
        ggtitle(NULL)
      
      p3 <- ggplot(ROC_QTL_archetecture_pooled, aes(False_positive_rate,Proportion_of_QTL_detected_weighted, color = factor(color, labels = labels))) + 
        #geom_point(size = 3, alpha = 0.8, shape=1) + 
        geom_line(size = 2, alpha = 0.8, linetype="dotted") + 
        #geom_point(data=ROC_QTL_archetecture_pooled[which(ROC_QTL_archetecture_pooled[,2]>0 | ROC_QTL_archetecture_pooled[,3]>0),], aes(False_positive_rate, Proportion_of_QTL_detected_weighted, color = factor(color, labels = labels)), size = 3, alpha = 0.8) + 
        geom_line(data=ROC_QTL_archetecture_pooled[which(ROC_QTL_archetecture_pooled[,2]>0 | ROC_QTL_archetecture_pooled[,3]>0),], aes(False_positive_rate, Proportion_of_QTL_detected_weighted, color = factor(color, labels = labels)), size = 2, alpha = 1) + 
        scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) + 
        xlab("False positive rate") +
        ylab("Proportion of QTL detected weighted by effect size") +
        guides(color=guide_legend(title="")) +
        theme(legend.position = c(0.65, 0.2)) +
        theme(axis.title = element_text(size = 30)) +
        theme(axis.text=element_text(size=30)) +
        theme(legend.text=element_text(size=30)) +
        theme(legend.background = element_rect(fill = "transparent")) +
        theme(text = element_text(size=30)) +
        ggtitle(NULL)
      
      p4 <- ggplot(ROC_QTL_archetecture_pooled, aes(False_positive_rate,Proportion_of_QTL_detected_weighted, color = factor(color, labels = labels))) + 
        #geom_point(size = 3, alpha = 0.8, shape=1) + 
        geom_line(size = 2, alpha = 0.8, linetype="dotted") + 
        #geom_point(data=ROC_QTL_archetecture_pooled[which(ROC_QTL_archetecture_pooled[,2]>0 | ROC_QTL_archetecture_pooled[,3]>0),], aes(False_positive_rate, Proportion_of_QTL_detected_weighted, color = factor(color, labels = labels)), size = 3, alpha = 0.8) + 
        geom_line(data=ROC_QTL_archetecture_pooled[which(ROC_QTL_archetecture_pooled[,2]>0 | ROC_QTL_archetecture_pooled[,3]>0),], aes(False_positive_rate, Proportion_of_QTL_detected_weighted, color = factor(color, labels = labels)), size = 2, alpha = 1) + 
        scale_x_continuous(limits = c(0,0.105), breaks = seq(0,0.1,0.02)) + 
        #scale_y_continuous(limits = c(0,1)) + 
        xlab("False positive rate") +
        ylab("Proportion of QTL detected weighted by effect size") +
        guides(color=guide_legend(title="")) +
        theme(legend.position = c(0.65, 0.2)) +
        theme(axis.title = element_text(size = 30)) +
        theme(axis.text=element_text(size=30)) +
        theme(legend.text=element_text(size=30)) +
        theme(legend.background = element_rect(fill = "transparent")) +
        theme(text = element_text(size=30)) +
        ggtitle(NULL)
      
      setwd(paste0(OutPath, "Figures/ROC/", title, "/"))
      
      png(paste0(title[1], "_", DirectionName, "_Original.png"), width = 1500, height = 760, units = "px", pointsize = 20)
      multiplot(p1, p3, layout = matrix(c(1,2), ncol=2))
      dev.off()
      
      png(paste0(title[1], "_", DirectionName, "_Zoomed.png"), width = 1500, height = 760, units = "px", pointsize = 20)
      multiplot(p2, p4, layout = matrix(c(1,2), ncol=2))
      dev.off()
      
      png(paste0(title[1], "_", DirectionName, "_Poster.png"), width = 800, height = 750, units = "px", pointsize = 20)
      multiplot(p4, layout = matrix(c(1), ncol=1))
      dev.off()
      
      rm(ROC_QTL_archetecture_pooled)
    }
  }
} else if (Type %in% c("Replication", "Direction")){
  if (Type == "Replication"){
    NsExpRep <- NsExpRep
  } else{
    NsExpRep <- paste0(1, Directions)
  }
  for (n_exp_rep in NsExpRep){
    for (sim_rep_id in 1:NSimRep){
      setwd(paste0(OutPath, "Simulations/", filename, "/SimRep", sim_rep_id, "/", "ExpRep", n_exp_rep, "/")) ## set directory to experimental replication
      ROC <- read.table("ROC_TwoTimepoints.txt", sep = ",",  header = T, stringsAsFactors = F)
      colnames(ROC)<- c("Proportion_of_genetic_variance_in_gen_1", 
                        "Proportion_of_QTL_detected_weighted", 
                        "False_positive_rate",
                        "Window_size",
                        "Computational_method",
                        "Numbers_of_timepoints_sampled")
      ROC <- ROC[which(ROC$WindowSize==WindowSize & ROC$ComputationalMethod==ComputationalMethod),]
      if(sim_rep_id==1){
        ROC_temp <- ROC[,1:3]
      } else {
        ROC_temp <- ((sim_rep_id-1)*ROC_temp + ROC[,1:3])/sim_rep_id
      }
    }
    color <- which(NsExpRep == n_exp_rep)
    ROC_temp <- cbind(ROC_temp, as.character(n_exp_rep), color)
    if(n_exp_rep==NsExpRep[1]){
      ROC_final <- ROC_temp
    } else {
      ROC_final <- rbind(ROC_final, ROC_temp)
    }
  }
  p1 <- ggplot(ROC_final, aes(False_positive_rate,Proportion_of_genetic_variance_in_gen_1, color = factor(color, labels = NsExpRep))) + 
    geom_point(size = 3, alpha = 0.8) + 
    xlab("False positive rate") +
    ylab("Proportion of genetic variance in the first generation") +
    ## scale_colour_manual("Methods", values = c("blue", "red", "green", "violet")) +
    theme(legend.position="none") +
    theme(axis.title = element_text(size = 30)) +
    theme(axis.text=element_text(size=30)) +
    theme(text = element_text(size=30)) +
    ggtitle(NULL)
  
  
  
  p2 <- ggplot(ROC_final, aes(False_positive_rate,Proportion_of_genetic_variance_in_gen_1, color = factor(color, labels = NsExpRep))) + 
    geom_point(size = 3, alpha = 0.8) + 
    xlim(0, 0.05) + 
    ylim(0, max(ROC_final$Proportion_of_genetic_variance_in_gen_1[which(ROC_final$False_positive_rate<=0.06)])) + 
    xlab("False positive rate") +
    ylab("Proportion of genetic variance in the first generation") +
    ## scale_colour_manual("Methods", values = c("blue", "red", "green", "violet")) +
    theme(legend.position="none") +
    theme(axis.title = element_text(size = 30)) +
    theme(axis.text=element_text(size=30)) +
    theme(text = element_text(size=30)) +
    ggtitle(NULL)
  
  p3 <- ggplot(ROC_final, aes(False_positive_rate,Proportion_of_QTL_detected_weighted, color = factor(color, labels = NsExpRep))) + 
    geom_point(size = 3, alpha = 0.8) + 
    xlab("False positive rate") +
    ylab("Proportion of QTL detected weighted by effect size") +
    guides(color=guide_legend(title=title)) +
    theme(legend.position = c(0.65, 0.2)) +
    theme(axis.title = element_text(size = 30)) +
    theme(axis.text=element_text(size=30)) +
    theme(legend.text=element_text(size=30)) +
    theme(legend.background = element_rect(fill = "transparent")) +
    theme(text = element_text(size=30)) +
    ggtitle(NULL)
  
  
  p4 <- ggplot(ROC_final, aes(False_positive_rate,Proportion_of_QTL_detected_weighted, color = factor(color, labels = NsExpRep))) + 
    geom_point(size = 3, alpha = 0.8) + 
    xlim(0, 0.05) + 
    ylim(0, max(ROC_final$Proportion_of_QTL_detected_weighted[which(ROC_final$False_positive_rate<=0.06)])) + 
    xlab("False positive rate") +
    ylab("Proportion of QTL detected weighted by effect size") +
    guides(color=guide_legend(title=title)) +
    theme(legend.position = c(0.65, 0.2)) +
    theme(axis.title = element_text(size = 30)) +
    theme(axis.text=element_text(size=30)) +
    theme(legend.text=element_text(size=30)) +
    theme(legend.background = element_rect(fill = "transparent")) +
    theme(text = element_text(size=30)) +
    ggtitle(NULL)
  
  setwd(paste0(OutPath, "Figures/ROC/", Type, "/"))
  
  png(paste0(title, "_Original.png"), width = 1500, height = 760, units = "px", pointsize = 20)
  multiplot(p1, p3, layout = matrix(c(1,2), ncol=2))
  dev.off()
  
  png(paste0(title, "_Zoomed.png"), width = 1500, height = 760, units = "px", pointsize = 20)
  multiplot(p2, p4, layout = matrix(c(1,2), ncol=2))
  dev.off()
  
} else if (Type %in% c("ComputationalMethod")){
  for (Direction in Directions){
    for (sim_rep_id in 1:NSimRep){
      setwd(paste0(OutPath, "Simulations/", filename, "/SimRep", sim_rep_id, "/", "ExpRep", Direction, NExpRep, "/")) ## set directory to experimental replication
      ROC <- read.table("ROC_TwoTimepoints.txt", sep = ",",  header = T, stringsAsFactors = F)
      temp <- read.table("ROC_AllTimepoints.txt", sep = ",",  header = T, stringsAsFactors = F)
      ROC <- rbind(ROC, temp)
      ROC[,5]<-paste0(ROC[,5], "_", ROC[,6])
      ROC <- ROC[,-6]
      colnames(ROC)<- c("Proportion_of_genetic_variance_in_gen_1", 
                        "Proportion_of_QTL_detected_weighted", 
                        "False_positive_rate",
                        "Window_size",
                        "Computational_method")
      if(sim_rep_id==1){
        ROC_temp <- ROC[,1:3]
      } else {
        ROC_temp <- ((sim_rep_id-1)*ROC_temp + ROC[,1:3])/sim_rep_id
      }
    }
    color <- sapply(ROC[,5], function(x) {which(unique(ROC[,5])==x)})
    labels <- c("D (2 Timepoints)", 
               "WFABC (2 Timepoints)", 
               "ApproxWF (2 Timepoints)",
               "WFABC (10 Timepoints)",
               "ApproxWF (10 Timepoints)")
    Computational_method <-ROC[,5]
    ROC_final <- cbind(ROC_temp, Computational_method, color)
    rm(Computational_method)
    p1 <- ggplot(ROC_final, aes(False_positive_rate,Proportion_of_genetic_variance_in_gen_1, color = factor(color, labels = labels))) + 
      geom_point(size = 5, alpha = 0.8) + 
      xlab("False positive rate") +
      ylab("Proportion of genetic variance in the first generation") +
      ## scale_colour_manual("Methods", values = c("blue", "red", "green", "violet")) +
      theme(legend.position="none") +
      theme(axis.title = element_text(size = 30)) +
      theme(axis.text=element_text(size=30)) +
      theme(text = element_text(size=30)) +
      ggtitle(NULL)
    
    p1 <- ggplot(ROC_final, aes(False_positive_rate, Proportion_of_genetic_variance_in_gen_1, color = factor(color, labels = labels))) + 
      #geom_point(size = 3, alpha = 0.8, shape=1) + 
      geom_line(size = 2, alpha = 0.8, linetype="dotted") + 
      #geom_point(data=ROC_final[which(ROC_final[,1]>0 | ROC_final[,3]>0),], aes(False_positive_rate, Proportion_of_genetic_variance_in_gen_1, color = factor(color, labels = labels)), size = 3, alpha = 0.8) + 
      geom_line(data=ROC_final[which(ROC_final[,1]>0 | ROC_final[,3]>0),], aes(False_positive_rate, Proportion_of_genetic_variance_in_gen_1, color = factor(color, labels = labels)), size = 2, alpha = 1) + 
      scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) + 
      xlab("False positive rate") +
      ylab("Proportion of genetic variance in the first generation") +
      theme(legend.position="none") +
      theme(axis.title = element_text(size = 30)) +
      theme(axis.text=element_text(size=30)) +
      theme(text = element_text(size=30)) +
      ggtitle(NULL)
    
    p2 <- ggplot(ROC_final, aes(False_positive_rate, Proportion_of_genetic_variance_in_gen_1, color = factor(color, labels = labels))) + 
      #geom_point(size = 3, alpha = 0.8, shape=1) + 
      geom_line(size = 2, alpha = 0.8, linetype="dotted") + 
      #geom_point(data=ROC_final[which(ROC_final[,1]>0 | ROC_final[,3]>0),], aes(False_positive_rate, Proportion_of_genetic_variance_in_gen_1, color = factor(color, labels = labels)), size = 3, alpha = 0.8) + 
      geom_line(data=ROC_final[which(ROC_final[,1]>0 | ROC_final[,3]>0),], aes(False_positive_rate, Proportion_of_genetic_variance_in_gen_1, color = factor(color, labels = labels)), size = 2, alpha = 1) + 
      scale_x_continuous(limits = c(0,0.105), breaks = seq(0,0.1,0.02)) + 
      scale_y_continuous(limits = c(0,1)) + 
      xlab("False positive rate") +
      ylab("Proportion of genetic variance in the first generation") +
      theme(legend.position="none") +
      theme(axis.title = element_text(size = 30)) +
      theme(axis.text=element_text(size=30)) +
      theme(text = element_text(size=30)) +
      ggtitle(NULL)
    
    p3 <- ggplot(ROC_final, aes(False_positive_rate,Proportion_of_QTL_detected_weighted, color = factor(color, labels = labels))) + 
      #geom_point(size = 3, alpha = 0.8, shape=1) + 
      geom_line(size = 2, alpha = 0.8, linetype="dotted") + 
      #geom_point(data=ROC_final[which(ROC_final[,2]>0 | ROC_final[,3]>0),], aes(False_positive_rate, Proportion_of_QTL_detected_weighted, color = factor(color, labels = labels)), size = 3, alpha = 0.8) + 
      geom_line(data=ROC_final[which(ROC_final[,2]>0 | ROC_final[,3]>0),], aes(False_positive_rate, Proportion_of_QTL_detected_weighted, color = factor(color, labels = labels)), size = 2, alpha = 1) + 
      scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) + 
      xlab("False positive rate") +
      ylab("Proportion of QTL detected weighted by effect size") +
      guides(color=guide_legend(title="")) +
      theme(legend.position = c(0.65, 0.2)) +
      theme(axis.title = element_text(size = 30)) +
      theme(axis.text=element_text(size=30)) +
      theme(legend.text=element_text(size=30)) +
      theme(legend.background = element_rect(fill = "transparent")) +
      theme(text = element_text(size=30)) +
      ggtitle(NULL)
    
    p4 <- ggplot(ROC_final, aes(False_positive_rate,Proportion_of_QTL_detected_weighted, color = factor(color, labels = labels))) + 
      #geom_point(size = 3, alpha = 0.8, shape=1) + 
      geom_line(size = 2, alpha = 0.8, linetype="dotted") + 
      #geom_point(data=ROC_final[which(ROC_final[,2]>0 | ROC_final[,3]>0),], aes(False_positive_rate, Proportion_of_QTL_detected_weighted, color = factor(color, labels = labels)), size = 3, alpha = 0.8) + 
      geom_line(data=ROC_final[which(ROC_final[,2]>0 | ROC_final[,3]>0),], aes(False_positive_rate, Proportion_of_QTL_detected_weighted, color = factor(color, labels = labels)), size = 2, alpha = 1) + 
      scale_x_continuous(limits = c(0,0.105), breaks = seq(0,0.1,0.02)) + 
      scale_y_continuous(limits = c(0,1)) + 
      xlab("False positive rate") +
      ylab("Proportion of QTL detected weighted by effect size") +
      guides(color=guide_legend(title="")) +
      theme(legend.position = c(0.65, 0.2)) +
      theme(axis.title = element_text(size = 30)) +
      theme(axis.text=element_text(size=30)) +
      theme(legend.text=element_text(size=30)) +
      theme(legend.background = element_rect(fill = "transparent")) +
      theme(text = element_text(size=30)) +
      ggtitle(NULL)
    
    setwd(paste0(OutPath, "Figures/ROC/", Type, "/"))
    
    png(paste0(title[1], "_", Direction, "_Original.png"), width = 1500, height = 760, units = "px", pointsize = 20)
    multiplot(p1, p3, layout = matrix(c(1,2), ncol=2))
    dev.off()
    
    png(paste0(title[1], "_", Direction, "_Zoomed.png"), width = 1500, height = 760, units = "px", pointsize = 20)
    multiplot(p2, p4, layout = matrix(c(1,2), ncol=2))
    dev.off()
    
    png(paste0(title[1], "_", Direction, "_Poster.png"), width = 800, height = 750, units = "px", pointsize = 20)
    multiplot(p4, layout = matrix(c(1), ncol=1))
    dev.off()
  }
}

