## This file is used to compile the slim output and also create input files for WFABC and ApproxWF.
## First define variables by selecting the experimental designs that you want to compare. This is necessary because compiled output from different experimental designs can generated from a single selection experiment simulation. (e.g. If you have two experimental replicates in the simulation, you can still include only one of them in power analysis.) Variable names are the same as in the shell and slim scripts before unless noted otherwise. 
## Run this script on server with the following command: nohup R CMD BATCH ./Compile.R  ./Compile.log > Compile.nohup &

rm(list=ls(all=TRUE))
library(data.table)

## Define variables
LCh=30000000 # length of chromosome; make sure that the simulation had LCh EQUAL to this
PopSize=1000 # make sure that the simulation had SampleSize EQUAL to this
SampleSize=50 # make sure that the simulation had SampleSize LARGER or EQUAL to this
NGen=10 # make sure that the simulation had NGen LARGER or EQUAL to this
NsExpRep=c(1 
             #, 2, 10
             ) # make sure that the simulation had NExpRep LARGER or EQUAL to each element in this
TimepointsSampled<-list(
  c(1:NGen),
  c(1, NGen)
) 
OutNames <- c(
  "AllTimepoints",
  "TwoTimepoints"
) # Just to differentiate output of all timepoints from two timepoints.
Directions<- c("Plus" 
               , "Minus"
) # Direction of selection; when doing both directions with repetition, the ROC table of single direction with repetition will not be outputted after the first replicate, so single direction with repetition should be reran individually if needed.
NSimRep=100 # make sure that the simulation had NSimRep LARGER or EQUAL to this
OutPath="/fs/cbsubscb10/storage/rl683/TemporalScan/Simulations/"
Filenames <- c("NQTL100"
               , "NQTL100_D0"
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
) # This is a list of directories under the "OutPath" directory that you have stored your selection simulation output at.


## Compile simulation outputs into clean dataframes for analyses ###################################################################
for (FilenameNumber in 22:length(Filenames)){
  print(FilenameNumber)
  filename <- Filenames[FilenameNumber]
  for (SimRepID in 1:NSimRep){
    for (NExpRep in 1:max(NsExpRep)){
      for (Direction in Directions){
        setwd(paste0(OutPath, filename, "/SimRep", SimRepID, "/", "ExpRep",  Direction, NExpRep, "/")) ## set directory to experimental replication
        for (i in 1:NGen){                                             # iterate through number of generations
          assign(paste0("Gen", i, "_Sample"), read.table(paste0("Gen", i, "_Sample.txt"), fill=T))
          n_SNP_per_gen <- which(get(paste0("Gen", i,"_Sample"))[,1]=="Genomes:") - 2 # get number of SNPs
          assign(paste0("Gen", i, "_Sample"), get(paste0("Gen", i, "_Sample"))[2:(n_SNP_per_gen+1),c(2,3,4,5,6,8,9)]) # trim off the idividual genotypes
          assign("temp", get(paste0("Gen", i, "_Sample"))) # in order to change colnames
          colnames(temp) <- c("ID_perm", paste0("type_",i), paste0("position_",i), paste0("es_",i), paste0("dominance_",i), paste0("origin_",i), paste0("freq_",i))
          if (i==1){
            assign("full_output_temp", temp)
          }else {
            assign("full_output_temp", merge(full_output_temp, temp, by="ID_perm", all=T))
          }  
        }
        n_SNP <- dim(full_output_temp)[1]
        full_output<-as.data.frame(matrix(NA, ncol=6+NGen,nrow=n_SNP)) # define a new full output for data cleaning
        colnames(full_output) <-c("ID_perm", "type",  "position", "es", "dominance", "origin", paste0("freq_",c(1:NGen)))
        full_output[,1]<-as.numeric(as.character(full_output_temp[,1])) # the first column is not duplicated
        my_max<-function(x){ # define this for the apply function
          return(max(x,na.rm=T)) 
        }
        
        ## Assign values to the new dataframe; the loop is necessary because there can be NAs in the first generation
        for (i in 2:6){ 
          full_output[,i]<-apply(full_output_temp[,seq(i, i+(NGen-1)*6, 6)], 1, my_max) 
        }
        full_output[,7:(6+NGen)]<-as.matrix(full_output_temp[,seq(7, 7+(NGen-1)*6, 6)])
        full_output[is.na(full_output)] <- 0 ## turn all NAs into 0
        
        ## Make all mutant alleles minor alleles (this is just for WFABC; can turn it off if WBABC is not used)
        for (i in 1:n_SNP){
          if (full_output[i,7] >SampleSize){
            full_output[i,7:(6+NGen)] <- SampleSize*2 - full_output[i,7:(6+NGen)]
            full_output[i,4] <- -full_output[i,4]
          }
        }
        full_output<-full_output[order(full_output[,4]),] ## order by effect size
        
        ## Make sure no mutation occurs at the same position; perhaps not necessary
        # while (anyDuplicated(full_output[,3]) > 0){
        #   full_output[duplicated(full_output[,3]), 3] <- full_output[duplicated(full_output[,3]), 3]+1  
        # }
        
        ## Generating Outputs
        for (j in 1:length(TimepointsSampled)){
          ## write the allele frequency output
          n_timepoints <- length(TimepointsSampled[[j]])
          freq_out <- full_output[,c(1:6,(6+TimepointsSampled[[j]]))]
          freq_out <- freq_out[apply(freq_out[,7:(6+n_timepoints)], 1, function(x) !all(x==0)),] # remove all SNPs that were never sampled
          n_SNP<-dim(freq_out)[1]
          write.table(freq_out, paste0("SampleFreq", OutNames[j],".txt"), sep=",", col.names=T, row.names = F, quote = F) # output the allele frequency
          
          ## write the WFABC input
          input_abc <- matrix(NA, ncol=n_timepoints, nrow = 2*n_SNP+2)
          input_abc[2, ] <- TimepointsSampled[[j]]
          for (i in 1:n_SNP){
            input_abc[2*i+1,] <- rep(SampleSize*2, n_timepoints)
            input_abc[2*i+2, ] <- as.matrix(freq_out[i,7:(6+n_timepoints)])
          }
          input_abc <- unname(input_abc)
          input_abc_new <- apply(input_abc, 1, paste, collapse = "," ) # for formatting
          input_abc_new[1] <- paste(n_SNP, n_timepoints, sep=" ")
          write.table(input_abc_new, paste0("ABCInput", OutNames[j],".txt"), sep=",", col.names=F, row.names = F, quote=FALSE, na = "")
          
          ## write the approxWF input
          input_app <- matrix(NA, nrow = n_timepoints+1, ncol = n_SNP+1)
          for (i in 1:n_timepoints){
            for (k in 1:n_SNP)
              input_app[i+1, k+1] <- paste0(input_abc[2*k+2, i], "/", input_abc[2*k+1, i])
          }
          input_app[1,-1]<-freq_out[,1]
          input_app[,1] <- c("time", TimepointsSampled[[j]])
          write.table(input_app,paste0("AppInput", OutNames[j],".txt"), quote = F, col.names=F, row.names = F)
        }
      }
    }
  }
}
