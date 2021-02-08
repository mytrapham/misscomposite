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
    - i: 1, 2 (this is used for parallel running of the do files on UCL's Myriad in order to create 2000 independent simulation repetitions for each simulation scenario);
3. ansim.do is a do file which analyses the simulation results.

### Results (output folder)
1. Datasets containing the simulation results are stored in files with naming format 'type'\_'case'\_'miss'\_'i'.dta; these files are analysed by ansim.do. 
2. Random-number states for the simulation runs of the complex composite endpoint are stored in files with naming format complex\_'case'\_'miss'\_rng_'i'.dta. These can be used to re-create simulation runs in which perfect prediction occurred during MI. 

### Log files (log folder)
Log files of the simulation study are stored with naming format 'type'\_'case'\_'miss'\_'i'.smcl

## Re-analysis of the TOPPS trial (misscomposite/topps_reanalysis)
The directory misscomposite/topps_reanalysis provides Stata code for the re-analysis of the TOPPS trials. Methods for handling missing values in the components of the composite endpoint are applied in the TOPPS data.

## References
<a id="1">[1]</a> 
Jann B (2005). 
MOREMATA: Stata module (Mata) to provide various functions.
Statistical Software Components S455001, Boston College Department of Economics, revised 06 Dec 2020.
