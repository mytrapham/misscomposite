//Outcome types: 1-simple; 2-complex
//Cases: 1-case I (noint); 2-case II (intctr); 3-case III (inttrtctr)
//Missingness mechanisms: 1-MCAR; 2-MAR1; 3-MAR2

*** COMPLEX NOINT MCAR ***
********************************************************************************

version 15

cd /lustre/scratch/scratch/zcakf79/Stata_output/

clear *
macro drop _all

cap log close
log using complex_noint_mcar_2, smcl replace

set rng mt64 
set rngstream 2001
set seed 8357

local R = 1000
local n = 2000
local m = 10
local c = 15

tempname postseed
postfile `postseed' str2000(s1 s2) str1100(s3) using "complex_noint_mcar_rng_2.dta", replace

tempname sim
postfile `sim' int(rep) str28(method) float(mcra mderiv b0 bx se0 sex df0 dfx errcode) using complex_noint_mcar_2, replace
	   			   
qui	{
	noi _dots 0, title("Simulation running...")
	timer on 1	

	forval r = 1/`R' 	{
		clear
		noi _dots `r' 0
		
		post `postseed' (substr(c(rngstate),1,2000)) (substr(c(rngstate),2001,2000)) (substr(c(rngstate),4001,.))

		//Outcome types: 1-simple; 2-complex
		//Cases: 1-case I (noint); 2-case II (intctr); 3-case III (inttrtctr)
		//Missingness mechanisms: 1-MCAR; 2-MAR1; 3-MAR2
		compositesim, outcometype(2) case(1) missingness(1) n(`n') m(`m') burnin(`c')
		
		local errcode = .
		
		post `sim' (`r') ("Full data") (r(mcra)) (r(mderiv)) (r(beta0_fd)) (r(betax_fd)) ///
									   (r(se0_fd)) (r(sex_fd)) (r(df0_fd)) (r(dfx_fd)) (`errcode')
								 
		post `sim' (`r') ("CRA") (r(mcra)) (r(mderiv)) (r(beta0_cra)) (r(betax_cra)) ///
							     (r(se0_cra)) (r(sex_cra)) (r(df0_cra)) (r(dfx_cra)) (`errcode')
		
		post `sim' (`r') ("Deriv") (r(mcra)) (r(mderiv)) (r(beta0_deriv)) (r(betax_deriv)) ///
								   (r(se0_deriv)) (r(sex_deriv)) (r(df0_deriv)) (r(dfx_deriv)) (`errcode')
		
		post `sim' (`r') ("MI-CRA") (r(mcra)) (r(mderiv)) (r(beta0_micra)) (r(betax_micra)) ///
								    (r(se0_micra)) (r(sex_micra)) (r(df0_micra)) (r(dfx_micra)) (`errcode')
		
		post `sim' (`r') ("MI-Deriv") (r(mcra)) (r(mderiv)) (r(beta0_mideriv)) (r(betax_mideriv)) ///
									  (r(se0_mideriv)) (r(sex_mideriv)) (r(df0_mideriv)) (r(dfx_mideriv)) (`errcode')
		
		post `sim' (`r') ("MIC-main") (r(mcra)) (r(mderiv)) (r(beta0_micmain)) (r(betax_micmain)) ///
										       (r(se0_micmain)) (r(sex_micmain)) (r(df0_micmain)) (r(dfx_micmain)) (`errcode')
											   
		post `sim' (`r') ("MIC-x") (r(mcra)) (r(mderiv)) (r(beta0_micx)) (r(betax_micx)) ///
													   (r(se0_micx)) (r(sex_micx)) (r(df0_micx)) (r(dfx_micx)) (r(by1_err)) 
													   
		post `sim' (`r') ("MIC-x-z1") (r(mcra)) (r(mderiv)) (r(beta0_micxz1)) (r(betax_micxz1)) ///
													      (r(se0_micxz1)) (r(sex_micxz1)) (r(df0_micxz1)) (r(dfx_micxz1)) (r(by2_err))	 				
	}
	 
	timer off 1	

}

timer list

postclose `sim'
postclose `postseed'

log close