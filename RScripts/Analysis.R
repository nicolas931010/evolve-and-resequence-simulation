## This script is used to calculate power and false positive rate at different thresholds.
## First define variables by selecting the experimental designs that you want to compare. This is necessary because compiled output from different experimental designs can generated from a single selection experiment simulation. (e.g. If you have two experimental replicates in the simulation, you can still include only one of them in power analysis.) Variable names are the same as in the shell and slim scripts before unless noted otherwise. 
## Run this script on server with the following command: nohup R CMD BATCH ./Analysis.R  ./Analysis.log > Analysis.nohup & 

rm(list=ls(all=TRUE))
library(data.table)

LChr=30000000 # length of chromosome; make sure that the simulation had LChr EQUAL to this
PopSize=1000 # make sure that the simulation had SampleSize EQUAL to this
SampleSize=50 # make sure that the simulation had SampleSize LARGER or EQUAL to this
NGen=10 # make sure that the simulation had NGen LARGER or EQUAL to this
NsExpRep=c(1 
             #, 2, 5, 10
             ) # make sure that the simulation had NExpRep LARGER or EQUAL to each element in this
TimepointsSampled<-list(
  #c(1:NGen),
  c(1, NGen)
  ) # Use c(1:NGen) only when comparing computational methods
OutNames <- c(
  #"AllTimepoints",
  "TwoTimepoints"
  ) # Use AllTimepoints only when comparing computational methods
Directions<- c("Plus" 
               , "Minus"
) # Direction of selection; when doing both Directions with repetition, the ROC table of single Direction with repetition will not be outputted, so single Direction with repetition should be ran if needed
WindowSizes <- c(0
                  #, 1000
                  #, 10000, 100000, 1000000
                 )
ComputationalMethods <- c("D"
                           #,"WFABC", "ApproxWF"
                          )
NSimRep=100 # make sure that the simulation had NSimRep LARGER or EQUAL to this
NROCPoints <- 100 # number of points to include in the ROC curve
OutPath="/fs/cbsubscb10/storage/rl683/TemporalScan/Simulations/" 
Filenames <- c("NQTL100"
               ,"NQTL100_D0"
               , "NQTL100_ESDistE"
               , "NQTL100_FreqLow5"
               , "NQTL10_D0"
               , "NQTL10_ESDistE"
               , "NQTL10_FreqLow5"
               , "NQTL10"
               , "NQTL100_Clustered"
               , "NQTL100_D1"
               , "NQTL100_FreqHigh5"
               , "NQTL10_Clustered"
               , "NQTL10_D1"
               , "NQTL10_FreqHigh5"
               , "NQTL10_EpiSce1"
               , "NQTL10_EpiSce2"
               , "NQTL10_EpiSce3"
               , "NQTL10_EpiSce4"
               , "NQTL10_EpiSce5"
               , "NQTL10_EpiSce6"
               , "NQTL10_FewerSNPs"
               , "NQTL10_D25"
               , "NQTL10_D75"
               , "NQTL100_D25"
               , "NQTL100_D75"
               , "NQTL200"
               , "NQTL2"
               , "NQTL20"
)

## Create ROC tables ###################################################################################################
########################################################################################################################
for (FilenameNumber in 22:length(Filenames)){ 
  Filename <- Filenames[FilenameNumber]
  for (SimRepID in 1:NSimRep){
    ## Get individual genotypes at true QTLs and individual trait values in the first generation; it's the same for all experimental replicates
    setwd(paste0(OutPath, Filename, "/SimRep", SimRepID, "/", "ExpRepPlus1/")) ## set directory to experimental replication
    data <- read.table("Gen1_Full.txt", fill = T, stringsAsFactors = F)
    gen_1_n_SNP <- which(data[,1]=="Genomes:") - 3
    gen_1_mutations <- data[2:(gen_1_n_SNP+1),c(1,2,3,4,5,6,8,9)]
    colnames(gen_1_mutations) <- c("ID_temp", "ID_perm", "type",  "position", "es", "dominance", "origin","freq")
    rm(data)
    gen_1_QTL_ID_temp <- as.numeric(gen_1_mutations[which(gen_1_mutations$es!=0),1])
    gen_1_QTL_ID_perm <- as.numeric(gen_1_mutations[which(gen_1_mutations$es!=0),2])
    gen_1_QTL_pos <- as.numeric(gen_1_mutations[which(gen_1_mutations$es!=0),4])
    gen_1_QTL_effect_size <- as.numeric(gen_1_mutations[which(gen_1_mutations$es!=0),5])
    gen_1_QTL_dominance <- as.numeric(gen_1_mutations[which(gen_1_mutations$es!=0),6])
    n_col<-max(count.fields("Gen1_Full.txt"))
    data <- read.table("Gen1_Full.txt", fill=T, col.names= 1:n_col,stringsAsFactors = F)
    gen_1_genomes <- data[(gen_1_n_SNP+4):(dim(data)[1]), ]
    gen_1_genomes <- unname(gen_1_genomes[,-(1:2)])
    gen_1_genomes <- data.matrix(gen_1_genomes)
    gen_1_QTL_phenotype <- matrix(0, ncol = length(gen_1_QTL_ID_perm), nrow = PopSize)
    colnames(gen_1_QTL_phenotype) <- gen_1_QTL_ID_perm
    for (i in 1:PopSize){
      for (j in 1:length(gen_1_QTL_ID_temp)){
        if (any(gen_1_genomes[2*i-1,]==gen_1_QTL_ID_temp[j], na.rm = T) & any(gen_1_genomes[2*i,]==gen_1_QTL_ID_temp[j], na.rm = T)){
          gen_1_QTL_phenotype[i,j]<-2*gen_1_QTL_effect_size[j]
        }
        else if (any(gen_1_genomes[2*i-1,]==gen_1_QTL_ID_temp[j], na.rm = T) | any(gen_1_genomes[2*i,]==gen_1_QTL_ID_temp[j], na.rm = T)){
          gen_1_QTL_phenotype[i,j]<-2*gen_1_QTL_effect_size[j]*gen_1_QTL_dominance[j]
        }
      }
    }
    gen_1_phenotype <- apply(gen_1_QTL_phenotype, 1, sum)
    
    ## Loop over different time sampling design (e.g. two points vs. full trajectory)
    for (TimepointsSampled_index in 1:length(TimepointsSampled)) {
      n_timepoints <- length(TimepointsSampled[[TimepointsSampled_index]])
      out_name <- OutNames[TimepointsSampled_index]
      if (exists("full_output")){ # remove "full_output" for the later step
        rm(full_output)
      } 
      ## Proceed to experimental replicates
      for (NExpRep in 1:max(NsExpRep)){
        for (Direction in Directions){
          setwd(paste0(OutPath, Filename, "/SimRep", SimRepID, "/", "ExpRep", Direction, NExpRep, "/")) ## set directory to experimental replication
          full_output_temp <-read.table(paste0("SampleFreq", out_name, ".txt"), sep = ",", header = T, stringsAsFactors = F)
          
          ## Contruct a dataframe consisting only of mutation info and D (the absolute value) and p-value
          ## then merge different exp_reps, using NExpRep to weigh the mean of D and p-value
          ## then we can use these D and p-values directly in later steps
          
          full_output_temp <- cbind(full_output_temp[,c(1:6)], full_output_temp[,6+n_timepoints]-full_output_temp[,7])
          if (Direction == "Plus"){
            full_output_temp[,7] <- full_output_temp[,7]/(2*SampleSize)
          } else {
            full_output_temp[,7] <- -full_output_temp[,7]/(2*SampleSize)
          }
          colnames(full_output_temp)[7]<-"D"
          full_output_temp <- cbind(full_output_temp, s_ABC=NA, s_App=NA)
          n_SNP <- dim(full_output_temp)[1]
          #########
          for (computational_method in ComputationalMethods){
            if (computational_method == "WFABC"){
              ## calculate p-value from WFABC output
              post_s <- read.table(paste0("ABCOutput", out_name, ".txt"))
              post_s <- as.data.frame(post_s)
              s_ABC<-rep(0,n_SNP) ## Attention to the n_SNP
              for (i in 1:n_SNP){
                s_ABC[i] <-mean(unlist(post_s[i,])) ## that's not the actual p-value
              }
              if (Direction == "Plus"){
                full_output_temp$s_ABC <- s_ABC 
              } else {
                full_output_temp$s_ABC <- -s_ABC
              }
            } else if (computational_method == "ApproxWF"){
              ## calculate p-value from ApproxWF output
              mcmc <- fread(paste0("AppOutput", out_name,".txt"))
              mcmc<-data.frame(mcmc)
              mcmc <- mcmc[10001:length(mcmc[,1]),]
              s_App <- rep(0,n_SNP) 
              for (i in 1:n_SNP){
                s_App[i] <- mean(mcmc[,i+2])  ## that's not the actual p-value ## Attention, might need to exclude first four columns depending on the ApproxWF setting!
              }
              if (Direction == "Plus"){
                full_output_temp$s_App <- s_App ## this is the proportion of posterior distribution larger than 0
              } else {
                full_output_temp$s_App <- -s_App
              }
            }
          }
          ## Merge different experimental replicates
          if (!exists("full_output")){
            full_output <- full_output_temp
          } else{
            full_output_temp <- merge(full_output, full_output_temp, by="ID_perm", all=T)
            n_SNP<-dim(full_output_temp)[1]
            full_output<-as.data.frame(matrix(NA, ncol=9,nrow=n_SNP)) # define a new full output for data cleaning
            full_output[,1]<-as.numeric(as.character(full_output_temp[,1])) # the first column is not duplicated
            my_max<-function(x){ # define this for the apply function
              return(max(x,na.rm=T)) 
            }
            ## Define weight to weigh the summary statistic when taking average
            weight <- length(Directions)*NExpRep + which(Directions==Direction) - length(Directions)
            ## Assign values to the new dataframe; the loop is necessary because there can be NAs in either dataframes
            for (i in 2:6){ 
              full_output[,i]<-apply(full_output_temp[,c(i, i+8)], 1, my_max) 
            }
            full_output_temp[,c(7, 15)][is.na(full_output_temp[,c(7, 15)])] <- 0 ## turn all NAs in D into 0
            full_output[,7]<-((weight-1)*full_output_temp[,7]+full_output_temp[,15])/weight # Calculate average D
            
            if("WFABC" %in% ComputationalMethods){
              full_output_temp[,c(8, 16)][is.na(full_output_temp[,c(8, 16)])] <- 1 ## turn all NAs in s_ABC into 1
              full_output[,8]<-((weight-1)*full_output_temp[,8]+full_output_temp[,16])/weight # Calculate average p-value
            }
            
            if("ApproxWF" %in% ComputationalMethods){
              full_output_temp[,c(9, 17)][is.na(full_output_temp[,c(9, 17)])] <- 1 ## turn all NAs in s_App into 1
              full_output[,9]<-((weight-1)*full_output_temp[,9]+full_output_temp[,17])/weight # Calculate average p-value
            }
            colnames(full_output) <-c("ID_perm", "type",  "position", "es", "dominance", "origin", "D", "s_ABC", "s_App")
          }
          if (NExpRep %in% NsExpRep){
            
            ## Define the sampled QTL positions, effect sizes, true interval
            sampled_SNP_effect_size <- full_output[,4]
            sampled_SNP_ID_perm <- full_output[,1]
            sampled_SNP_pos <- full_output[,3]
            sampled_QTL_pos <- sampled_SNP_pos[which(sampled_SNP_effect_size!=0)]
            sampled_QTL_ID <- sampled_SNP_ID_perm[which(sampled_SNP_effect_size!=0)]
            sampled_QTL_effect_size <- sampled_SNP_effect_size[which(sampled_SNP_effect_size!=0)]
            
            ######################################################################################################################
            
            ## Quantify selection by SNP
            for (computational_method in ComputationalMethods){
              if (out_name=="AllTimepoints" & computational_method == "D"){
              } else {
                if (computational_method == "D"){
                  ## get abs_D from full_output
                  D <- full_output$D
                  abs_D <- abs(D)
                } else if (computational_method == "WFABC"){
                  ## get WFABC p-value from from full_output
                  abs_s_ABC<-abs(full_output$s_ABC)
                } else if (computational_method == "ApproxWF"){
                  ## get log ApproxWF p-value from full_output
                  abs_s_App <- abs(full_output$s_App)
                }
                
                ## Loop over different window sizes for window-based analyses
                for (window_size in WindowSizes){ 
                  ## First define true intervals and SNPs in true intervals
                  gen_1_SNP_ID_perm <- as.numeric(gen_1_mutations[,2])
                  gen_1_SNP_pos <- as.numeric(gen_1_mutations[,4])
                  true_interval <- cbind(gen_1_QTL_pos-window_size/2, gen_1_QTL_pos+window_size/2)
                  true_interval[which(true_interval[,1] < 0), 1] <- 0
                  true_interval[which(true_interval[,2] > (LChr-1)), 2] <- (LChr-1)
                  gen_1_SNP_in_true_interval_pos <- gen_1_SNP_pos[inrange(gen_1_SNP_pos, true_interval[,1], true_interval[,2])]
                  
                  ## Create empty vectors to contain power and false positive rate
                  Power_gen_1 <- rep(0,NROCPoints)
                  Power_percentage_QTL_detected <- rep(0, NROCPoints)
                  False_positive <- rep(0, NROCPoints)
                  
                  ## Calculate power and false positive rates at different threshold
                  for (i in 1:NROCPoints){
                    if (computational_method == "D"){
                      cutoff <- 1*(i-1)/NROCPoints
                      detected_pos <- sampled_SNP_pos[which(abs_D>=cutoff)]
                    } else if (computational_method == "WFABC"){
                      cutoff <- 1*(i-1)/NROCPoints
                      detected_pos <- sampled_SNP_pos[which(abs_s_ABC>=cutoff)]
                    } else if (computational_method == "ApproxWF"){
                      cutoff <- 1*(i-1)/NROCPoints
                      detected_pos <- sampled_SNP_pos[which(abs_s_App>=cutoff)]
                    }
                    detected_interval<-cbind(detected_pos-window_size/2,detected_pos+window_size/2)
                    detected_interval[which(detected_interval[,1] < 0), 1] <- 0
                    detected_interval[which(detected_interval[,2] > (LChr-1)), 2] <- (LChr-1)
                    detected_QTL_pos <- sampled_QTL_pos[inrange(sampled_QTL_pos, detected_interval[,1], detected_interval[,2])]
                    detected_QTL_ID_perm <- sampled_QTL_ID[inrange(sampled_QTL_pos, detected_interval[,1], detected_interval[,2])]
                    SNP_in_detected_interval_pos <- sampled_SNP_pos[inrange(sampled_SNP_pos, detected_interval[,1], detected_interval[,2])]
                    SNP_in_detected_interval_ID <- sampled_SNP_ID_perm[inrange(sampled_SNP_pos, detected_interval[,1], detected_interval[,2])]
                    
                    ## Genetic variance explained in the first generation
                    if (length(detected_QTL_ID_perm)==0){
                      Power_gen_1[i] <- 0
                      Power_percentage_QTL_detected[i] <- 0
                    } else {
                      detected_gen_1_QTL_phenotype <- data.frame(gen_1_QTL_phenotype[,which(colnames(gen_1_QTL_phenotype) %in% detected_QTL_ID_perm)])
                      detected_gen_1_phenotype <- apply(detected_gen_1_QTL_phenotype, 1, function(x) sum(x, na.rm = T))
                      Power_gen_1[i] <- cor(gen_1_phenotype, detected_gen_1_phenotype)^2
                      ## Percentage QTL detected weighted by effect size
                      Power_percentage_QTL_detected[i] <- sum(abs(sampled_QTL_effect_size[which(sampled_QTL_pos %in% detected_QTL_pos)]))/sum(abs(gen_1_QTL_effect_size)) 
                    }
                    ## False positive
                    False_positive[i] <-  (length(SNP_in_detected_interval_pos)-length(which(SNP_in_detected_interval_pos %in% gen_1_SNP_in_true_interval_pos)))/(length(sampled_SNP_pos)-length(which(sampled_SNP_pos %in% gen_1_SNP_in_true_interval_pos)))
                    
                  }
                  Power_gen_1[is.na(Power_gen_1)]<-0
                  Window_size <- rep(window_size, NROCPoints)
                  ROC_temp<-data.frame(Power_gen_1,Power_percentage_QTL_detected, False_positive, Window_size, computational_method, out_name)
                  colnames(ROC_temp)<- c("Proportion_of_genetic_variance_in_gen_1", 
                                         "Proportion_of_QTL_detected_weighted", 
                                         "False_positive_rate",
                                         "Window_size",
                                         "Computational_method",
                                         "Numbers_of_timepoints_sampled")
                  if (!exists("ROC")){
                    ROC<-ROC_temp
                  } else {
                    ROC<-rbind(ROC, ROC_temp)
                  }
                }
              }
            }
            ## Also, when doing selection in both Directions, the output only makes sense when results from both Directions are incorporated in it
            if (exists("ROC") & Direction==Directions[length(Directions)]){ 
              write.table(ROC, paste0("ROC_", out_name, ".txt"), sep=",", col.names=T, row.names = F, quote = F)
            } else if (exists("ROC") & NExpRep==1){
              write.table(ROC, paste0("ROC_", out_name, ".txt"), sep=",", col.names=T, row.names = F, quote = F)
            }
            ## The following is necessary because sometimes ROC is not created
            if (exists("ROC")){
              rm(ROC)
            }
          }
        }
      }
    }
  }
}
