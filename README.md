# A comparison of methods for analysing a binary composite endpoint with partially observed components in randomised controlled trials

## Simulation study (misscomposite/simulation)
The directory misscomposite/simulation provides Stata code and results for the simulation study which explores a set of methods for handling missing values in the three components of a simple/complex composite endpoint.

### Stata code (do folder)
There are 3 sets of files:
1. composite.ado is a programme which executes a single simulation run (installation of moremata [[1]](#1) is required);
2. Files with naming format 'type'\_'case'\_'miss'\_'i'.do are do files which create 1000 repetitions of the simulation programme in composite.ado;
    -  type: simple, complex;
    - case: noint (case I), intctr (case II), inttrtctr (case III);
    - miss: mcar, mar1, mar2;
    - i: 1, 2 (this is used for parallel running of the do files on UCL's Myriad High Performance Computing cluster in order to create 2000 independent simulation repetitions for each simulation scenario);
3. ansim.do is a do file which analyses the simulation results (installation of simsum [[2]](#2) is required).

### Results (output folder)
1. Datasets containing the simulation results are stored in files with naming format 'type'\_'case'\_'miss'\_'i'.dta; these files are analysed by ansim.do. 
2. Random-number states for the simulation runs of the complex composite endpoint are stored in files with naming format complex\_'case'\_'miss'\_rng_'i'.dta. These can be used to re-create simulation runs in which perfect prediction occurred during MI. 

### Log files (log folder)
Log files of the simulation study are stored with naming format 'type'\_'case'\_'miss'\_'i'.smcl

## Re-analysis of the TOPPS trial (misscomposite/topps_reanalysis)
The directory misscomposite/topps_reanalysis provides Stata code (topps_reanalysis.do) for the re-analysis of the TOPPS trials (installation of ice [[3]](#3) is required). Methods for handling missing values in the components of the composite endpoint are applied to the TOPPS data. 

## References
<a id="1">[1]</a> 
Jann B (2005). 
MOREMATA: Stata module (Mata) to provide various functions.
Statistical Software Components S455001, Boston College Department of Economics, revised 06 Dec 2020.

<a id="2">[2]</a> 
White IR (2010). 
simsum: Analyses of simulation studies including Monte Carlo error.
The Stata Journal, 10(3):369–385.

<a id="3">[3]</a> 
Royston P, White IR (2011). 
Multiple imputation by chained equations (MICE): Implementation in Stata
Journal of Statistical Software, 45(4), 2011.
