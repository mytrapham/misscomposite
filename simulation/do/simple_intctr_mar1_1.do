//Outcome types: 1-simple; 2-complex
//Cases: 1-case I (noint); 2-case II (intctr); 3-case III (inttrtctr)
//Missingness mechanisms: 1-MCAR; 2-MAR1; 3-MAR2

*** SIMPLE INTCTR MAR1 ***
********************************************************************************

version 15

cd /lustre/scratch/scratch/zcakf79/Stata_output/

clear *
macro drop _all

cap log close
log using simple_intctr_mar1_1, smcl replace

set rng mt64 
set rngstream 1
set seed 4499

local R = 1000
local n = 2000
local m = 10
local c = 15

tempname sim
postfile `sim' int(rep) str28(method) float(mcra mderiv b0 bx se0 sex df0 dfx) using simple_intctr_mar1_1, replace
	   			   
qui	{
	noi _dots 0, title("Simulation running...")
	timer on 1	

	forval r = 1/`R' 	{
		clear
		noi _dots `r' 0	
		
		//Outcome types: 1-simple; 2-complex
		//Cases: 1-case I (noint); 2-case II (intctr); 3-case III (inttrtctr)
		//Missingness mechanisms: 1-MCAR; 2-MAR1; 3-MAR2
		compositesim, outcometype(1) case(2) missingness(2) n(`n') m(`m') burnin(`c')
								   
		post `sim' (`r') ("Full data") (r(mcra)) (r(mderiv)) (r(beta0_fd)) (r(betax_fd)) ///
									   (r(se0_fd)) (r(sex_fd)) (r(df0_fd)) (r(dfx_fd))
								 
		post `sim' (`r') ("CRA") (r(mcra)) (r(mderiv)) (r(beta0_cra)) (r(betax_cra)) ///
							     (r(se0_cra)) (r(sex_cra)) (r(df0_cra)) (r(dfx_cra))
		
		post `sim' (`r') ("Deriv") (r(mcra)) (r(mderiv)) (r(beta0_deriv)) (r(betax_deriv)) ///
								   (r(se0_deriv)) (r(sex_deriv)) (r(df0_deriv)) (r(dfx_deriv))
		
		post `sim' (`r') ("MI-CRA") (r(mcra)) (r(mderiv)) (r(beta0_micra)) (r(betax_micra)) ///
								    (r(se0_micra)) (r(sex_micra)) (r(df0_micra)) (r(dfx_micra))
		
		post `sim' (`r') ("MI-Deriv") (r(mcra)) (r(mderiv)) (r(beta0_mideriv)) (r(betax_mideriv)) ///
									  (r(se0_mideriv)) (r(sex_mideriv)) (r(df0_mideriv)) (r(dfx_mideriv))
		
		post `sim' (`r') ("MIC-main") (r(mcra)) (r(mderiv)) (r(beta0_micmain)) (r(betax_micmain)) ///
										       (r(se0_micmain)) (r(sex_micmain)) (r(df0_micmain)) (r(dfx_micmain))
											   
		post `sim' (`r') ("MIC-x") (r(mcra)) (r(mderiv)) (r(beta0_micx)) (r(betax_micx)) ///
													   (r(se0_micx)) (r(sex_micx)) (r(df0_micx)) (r(dfx_micx))	
													   
		post `sim' (`r') ("MIC-x-z1") (r(mcra)) (r(mderiv)) (r(beta0_micxz1)) (r(betax_micxz1)) ///
													      (r(se0_micxz1)) (r(sex_micxz1)) (r(df0_micxz1)) (r(dfx_micxz1))	
	}
	 
	timer off 1	

}

timer list

postclose `sim'

log close