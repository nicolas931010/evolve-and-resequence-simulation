# evolve-and-resequence-simulation
This is a customizable simulation tool for the power analysis of evolve and resequence experiments based on SLiM2.

IMPORTANT: This simulation tool is still in development. A few functionalities still need to be incorporated, most importantly the ability to change the heritability of the trait (which is set to 1 currently). Except for that, the script should be generally complete, especially for step 1-3. Step 4 and 5 might still be buggy but we are working on them. We also recommend to use your own R scripts for step 5. Please contact us at rl683@cornell.edu with any problems or questions. 

Note: There are some variables that will not be explicitly incorporated in this tool but should be fairly easy to implement if you dig into the slim script. These include selection modes other than truncating selection, population structure, and demographic history prior to the selection experiment. 


Here is our suggested workflow. Each step in the workflow has to be completed before the next can start
1. Create neutal populations representing the population prior to the selection experiment. To do this, read and edit ShellScripts/Burnin.sh , as well as SlimScripts/Burnin.slim ; then run Burnin.sh
2. Establish trait architecture and simulate selection experiment. For this, read ShellScripts/Selection.sh and SlimScripts/Selection.slim ; copy and rename Selection.sh and edit accordingly; then run the new script. This step needs to be repeated to simulate different combinations of trait architecture and selection experiment. They can be created using the same burnin. 
3. Compile SLiM outputs from each simulations in step 2, and (optionally) create input files for WFABC and AppWF. To do this, read, edit, and run RScripts/Compile.R 
4. Calculate power and false positive rate for each simulations. For this, read, edit, and run RScripts/Analysis.R
5. Plot ROC curve comparisons by combining the ROC tables from simulations of different trait architecture or experimental design given by step 4. You can use RScripts/PlotROC.R to do this. Alternatively, you can also create your own R scripts, since this step highly depends on what your specific comparison is. 


If using WFABC and/or ApproxWF. They need to be run between step 3 and 4. The output should be stored in the same directory as the input files. 
Various miscellaneous plots can be created using RScripts/PlotMisc.R and can help your analyses. More of these plots will be added to the script. 
