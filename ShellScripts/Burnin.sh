## This script is used to generate neutral populations used for the selection experiment.
## IMPORTANT: Change the output directory in the first step and the variables in the second step before running.
## Run the script on server using: nohup bash Burnin.sh > Burnin.nohup &
## Note: the number of generations to be included in the burnin process depends on the population size (typically ten times the population size) and needs to be changed in the Burnin.slim file. Please refer to the notes in the .slim file. 

## Step 1: Ceate directories to store the outputs. Make sure to change the directory name in the first line.

cd /fs/cbsubscb10/storage/rl683/TemporalScan/Simulations # Change this to a directory where you want to store all you simulation outputs.
mkdir Burnin
cd Burnin
for k in {1..100} # Number of simulation replicates that you want to create.
do 
    mkdir 'SimRep'$k
done
cd ..

## Step 2: Run burnin using SLiM 2. Variables inside the loop are all customizable and can be changed as desired. 

for k in {1..100} # Set the number of simulation replicates that you want to create.
do 
    ~/Program/SLiM/bin/slim \ # slim directory
    -d SimRepID=$k  \  # SimuRepID = Simulation Replicate ID
    -d Mu=2e-8 \ # Mu = mutation rate
    -d RecRate=1e-8 \ # RecRate = recombination rate (change the slim script if simulating multiple chromosomes)
    -d LCh=30000000 \ # LCh = length of chromosome 
    -d BurninSize=1000 \ # BurninSize = size of the burnin populations (IF THIS NEEDS TO BE CHANGED, MAKE SURE TO CHANGE MUTATION RATE, RECOMBINATION RATE, AND NUMBER OF GENERATIONS ACCORDINGLY) 
    -d "BurninPath='/fs/cbsubscb10/storage/rl683/TemporalScan/Simulations/Burnin/'" \ # BurninPath = path to the burnin files
    -d "BurninFilename='Burnin.txt'" \ # BurninFilename = name of the burnin file
    /fs/cbsubscb10/storage/rl683/TemporalScan/SlimScripts/Burnin.slim # Directory to the Burnin.slim file included in the simulation tool.
done
