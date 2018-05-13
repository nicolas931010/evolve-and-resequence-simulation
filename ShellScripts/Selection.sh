## This script is used to establish quantitative trait architectures and perform artificial selection experiments on the neutral populations generated using the Burnin.sh script.
## IMPORTANT: Change the output directory and in the first step and the variables in the second step before running.
## Run the script on server using: nohup bash Selection.sh > Selection.nohup &

## Step 1: Ceate directories to store the outputs. Make sure to change the directory names in the first three lines.

cd /fs/cbsubscb10/storage/rl683/TemporalScan/Simulations/ # Change this to a directory where you want to store all you simulation outputs.
mkdir NQTL10 # Change this to what you want to name this particular quantitative trait architectures and/or the experimental design that you are simulating.
cd NQTL10 # Same as above
for k in {1..100} # Number of simulation replicates that you want to create.
do
    mkdir 'SimRep'$k
    cd 'SimRep'$k
    for j in {1..10}
    do
        mkdir 'ExpRepPlus'$j
        mkdir 'ExpRepMinus'$j
    done
    cd ..
done
cd ..

## Step 2: Run burnin using SLiM 2. Variables inside the loop are all customizable and can be changed as desired. 
## IMPORTANT: Number of generations has to be changed in the slim script.
## It is recommended if you copy this shell script for each of the trait architecture x experimental combinations that you want to test and make changes in the new script. 
## For example, when simulating the scenario with 100 QTLs, copy this file and rename it as NQTL100.sh. Set the NQTL variable to be 100 and run it with "nohup bash NQTL100.sh > NQTL100.nohup &" 
## Note: The number of generations in the selection experiment can only be changed in the Selection.slim script.
## Depending on the population size, more or fewer burnin generations might have been needed. If that's the case, edit the Selection.slim file as instructed in the file. 

for k in {1..100} # Set the number of simulation replicates that you want to create.
do
    echo $k
    for j in {1..1} # Set the number of experimental replications. 
    do
        echo $j
        for i in {T,F} # Set the direction of selection (F if selecting the larger phenotype, T otherwise, T,F is both directions are selected)
        do
            echo $i
            ~/Program/SLiM/bin/slim \ # slim directory
            -d SimRepID=$k  \ # SimuRepID = Simulation Replicate ID
            -d ExpRepID=$j \ # ExpRepID = Experiment replicate ID
            -d Direction=$i \ # Direction = direction of selection; set this at the line "for i in {T,F}" above
            -d "BurninPath='/fs/cbsubscb10/storage/rl683/TemporalScan/Simulations/Burnin/'" \ # BurninPath = path to the burnin files
            -d "BurninFilename='Burnin.txt'" \ # BurninFilename = name of the burnin file
            -d LCh=30000000 \ # LCh = length of chromosome (CHANGE THIS ONLY WHEN THE BURNIN USES A DIFFERENT CHROMOSOME SIZE)
            -d RecRate=1e-8 \ # RecRate = recombination rate (change the slim script if simulating multiple chromosomes)
            -d SampleSize=50 \ # SampleSize = number of individuals to sample each generation
            -d NQTL=10 \ # NQTL = number of QTLs, (even number is recommended when the number is small)
            -d ESMean=1.0 \ # ESMean = absolute value of mean effect size
            -d "ESDist='f'" \ # ESDist = effect size distribution("f" for fixed or "e" for exponential), 
            -d LowFreq=F \ # LowFreq = starting frequency preference (T if selecting for lower frequency, F if selecting for higher frequency or random) 
            -d FreqBound=0.0 \ # FreqBound = frequency bound (0.0~0.5, 0.0 if starting frequency is random)
            -d LowerPosBound=0 \ # LowerPosBound = lower position bound (0 if random)
            -d UpperPosBound=29999999 \ # UpperPosBound = upper position bound (LCh-1 if random)
            -d D=0.5 \ # D = dominance coefficient (0.0~1.0, 1.0 for mutant being completely dominant and 0.0 for wildtype to being completely dominant.)
            -d Epistasis=F \ # Epistasis = F if there is no epistasis, T otherwise
            -d "EpiSce=c(0,1,2,1,2,3,2,3,4)" \ # EpiSce = epistasis scenario; vector of size 9; first element must be 0; needs to be defined if Epistasis ==T; the nine value corresponds to phenotypes of genotypes in the following order: c(aabb, Aabb, AAbb, aaBb, AaBb, AABb, aaBB, AaBB, AABB); e.g. c(0,1,2,1,2,3,2,3,4) when there is no epistasis.
            -d PopSize=1000 \ # PopSize = population size (can only downsample)
            -d SelectedSize=100 \ # SelectedSize = number of selected individuals in each generation
            -d "OutPath='/fs/cbsubscb10/storage/rl683/TemporalScan/Simulations/NQTL10/'" \ # OutPath = output path
            /fs/cbsubscb10/storage/rl683/TemporalScan/SlimScripts/Selection.slim # Directory to the Selection.slim file included in the simulation tool.
        done
    done
done
