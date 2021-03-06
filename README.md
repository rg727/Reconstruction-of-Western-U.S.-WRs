# A-Multi-Objective-Paleo-Informed-Reconstruction-of-Western-U.S.-Weather-Regimes-

Repository for Gupta, Steinschneider, & Reed, 2021, "A Multi-Objective Paleo-Informed Reconstruction of Western U.S. Weather Regimes Over the Past 600 Years" (in preparation)

This repository contains all code needed to run the optimization and generate the reconstruction. This experiment was run on THCUBE Cluster at Cornell University. If you have questions about running the code please contact rg727@cornell.edu.  

Note: One missing file is:

1) borg.R: The R wrapper for the Borg Multi-Objective Evolutionary Algorithm (http://borgmoea.org/), a licensed software. Please request access from the website. Alternatively, the optimization and reconstruction can be performed without Borg, using any standard multi-objective algorithm, though the same quality of results and convergence cannot be guaranteed for other algorithms. 

## Contents ##
* Code: Contains all code necessary to perform the optimization and reconstruction
* Data: Contains all code necessary to run the optimization

## Steps to Run for Performing the Optimization and Reconstruction (./code) ## 

The main file that performs the optimization and reconstruction is WR_Main.R. The main file references three reelvant scripts: "DirichReg_optim_reference.R" which performs the optimization routine, "DirichReg_predict_cv.R" which does the final cross-validated prediction for the chosen solution, and "DirichReg_reconstruct.R" which performs the reconstruction. 

