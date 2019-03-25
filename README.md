GM_DEBPCoD
==========
A life history model for _Globicephala melas_ (pilot whales) based on Dynamic Energy Budgets
----------

Code from: 
---------
Hin, V, Harwood, J. & De Roos, A.M.
_Bio-energetic modeling of medium-sized cetaceans shows high sensitivity to disturbance in seasons of low resource supply_
Ecological Applications 2019

Execution of the code requires the installation of EBTtool software package See: (https://staff.fnwi.uva.nl/a.m.deroos/EBT/Software/index.html) for code & installation

Code was executed on 2017 MacBook Pro running macOS High Sierra (10.13.6) with R version 3.5.1. Note that successful execution of the code might depend on your specific installation.

Code is released under GNU General Public License v3.0. A copy of the license is attached.

Please cite the above paper when reusing any part of the model or the code in a publication

Last modified: VH - 25 March, 2019

How to use
----------

* Obtain and install the EBTtool software package from https://staff.fnwi.uva.nl/a.m.deroos/EBT/Software/index.html. Make sure to define the appropriate environmental variables (EBTPATH and PATH) to run the program from the command line (terminal).
* In the terminal set the working directory to the folder that contains the model files 
* Invoke `make GM_DEBPCoD` to build the model
* Invoke `./GM_DEBPCoD single` to run the model with default parameter values as specified by `single.cvf`
* This creates 5 output files:
  1. For every daily timestep, `single.out` contains all the (environmental) output variables as defined by the routine DefineOutput() in `GM_DEBPCoD.c`
  2. `single.rep` is a report file that summarizes the parameter values and details of the integration of the executed run
  3. `single.esf` contains the values of the vectors that hold individual state variables (`i_state()`) and individual state constants (`i_const()`) for the last time step of the integration (which equals the end of life of the female)
  4. `single_R0.dat` reports life history statistics from the total lifetime of the female as specified by the lines 505-510 in `GM_DEBPCoD.c`
  5. `single.csb` normally hold the population state-output, but is empty since no population state output was recorded for this particular integration
  
* Run the R-file GM_DEBPCoD.R to import and plot the model output
* Source Figure_1.R to create the figure 1 from the manuscript
* Figures 2 & 3 requires changing the appropriate parameters in the parameter definition file (`single.cvf`)
* run `make allclean` to clean all model output files and c object files and executables
