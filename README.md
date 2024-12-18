This R code can be used to reproduce the simulations and analysis found in "Bootstrap Aggregated Designs for Generalized Linear Models" by Rios & Stufken (2024).
The following files are included:

1. BAGDesignRobustSim.R. This is the main file that includes all R functions required to find BAG designs. 
Lines 512-528 are used to specify the simulation parameters (see Section 4.1 of the paper). Once the parameters are specified, the simulation can be executed. The simulation was executed on GMU's HOPPER cluster.
If this is run on a single machine, it is recommended to use fewer cores (default is set to 30) and smaller values of n_sim (number of simulated datasets) and n_iter (number of iterations per simulated dataset).

2. SimResults.zip. This zip file contains csv files that store the output of BAGDesignRobustSim.R for all simulation parameters specified in Section 4.1. These results are passed to makerobustboxplots.R to generate
Figures 2 to 5.

3. makerobustboxplots.R. This R file takes csv files from SimResults.zip and uses them to generate Figures 2 to 5. To do this, unzip SimResults.zip and change line 1 of the R code to the directory where the csv files are
stored on your device.

4. BAGweightsonly.zip. This zip file contains the average BAG weights obtained for each of the 30 simulated datasets in Section 4. It was used to create Figure 6.

5. make_boxplots_weightsonly.R. This R file uses the csv files in BAGweightsonly.zip to create Figure 6.

6. AppendixBSimResults.zip. This zip file contains csv files that store the output of BAGDesignRobustSim.R for all of the simulations in Appendix B. These results may be passed to makerobustboxplots.R to generate all figures in Appendix B.

Our code requires the use of the R package matlib. If you encounter issues installing matlib on a Mac, please try and install XQuartz first, as per the following instructions from the package authors: “The functions that draw 3D graphs use the rgl package. On macOS, the rgl package requires that XQuartz be installed. After installing XQuartz, it’s necessary either to log out of and back into your macOS account or to reboot your Mac.”
