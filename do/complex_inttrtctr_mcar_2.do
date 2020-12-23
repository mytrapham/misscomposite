//Outcome type: 1-simple; 2-complex
//Case: 1-case I (noint); 2-case II (intctr); 3-case III (inttrtctr)
//Missingness: 1-MCAR; 2-MAR1; 3-MAR2

*** COMPLEX INTTRTCTR MCAR ***
********************************************************************************

version 15

cd /lustre/scratch/scratch/zcakf79/Stata_output/

clear *
macro drop _all

cap log close
log using complex_inttrtctr_mcar_2, smcl replace

set rng mt64 
set rngstream 2001
set seed 5515

local R = 1000
local n = 2000
local m = 10
local c = 15

tempname postseed
postfile `postseed' str2000(s1 s2) str1100(s3) using "complex_inttrtctr_mcar_rng_2.dta", replace

tempname sim
postfile `sim' int(rep) str28(method) float(mcra mderiv b0 btrt se0 setrt df0 dftrt errcode) using complex_inttrtctr_mcar_2, replace
	   			   
qui	{
	noi _dots 0, title("Simulation running...")
	timer on 1	

	forval r = 1/`R' 	{
		clear
		noi _dots `r' 0
		
		post `postseed' (substr(c(rngstate),1,2000)) (substr(c(rngstate),2001,2000)) (substr(c(rngstate),4001,.))
		
		//Outcome type: 1-simple; 2-complex
		//Case: 1-case I (noint); 2-case II (intctr); 3-case III (inttrtctr)
		//Missingness: 1-MCAR; 2-MAR1; 3-MAR2
		compositesim, outcometype(2) case(3) missingness(1) n(`n') m(`m') burnin(`c')
								   
		local errcode = .
		
		post `sim' (`r') ("Full data") (r(mcra)) (r(mderiv)) (r(beta0_fd)) (r(betatrt_fd)) ///
									   (r(se0_fd)) (r(setrt_fd)) (r(df0_fd)) (r(dftrt_fd)) (`errcode')
								 
		post `sim' (`r') ("CRA") (r(mcra)) (r(mderiv)) (r(beta0_cra)) (r(betatrt_cra)) ///
							     (r(se0_cra)) (r(setrt_cra)) (r(df0_cra)) (r(dftrt_cra)) (`errcode')
		
		post `sim' (`r') ("Deriv") (r(mcra)) (r(mderiv)) (r(beta0_deriv)) (r(betatrt_deriv)) ///
								   (r(se0_deriv)) (r(setrt_deriv)) (r(df0_deriv)) (r(dftrt_deriv)) (`errcode')
		
		post `sim' (`r') ("MI-CRA") (r(mcra)) (r(mderiv)) (r(beta0_micra)) (r(betatrt_micra)) ///
								    (r(se0_micra)) (r(setrt_micra)) (r(df0_micra)) (r(dftrt_micra)) (`errcode')
		
		post `sim' (`r') ("MI-Deriv") (r(mcra)) (r(mderiv)) (r(beta0_mideriv)) (r(betatrt_mideriv)) ///
									  (r(se0_mideriv)) (r(setrt_mideriv)) (r(df0_mideriv)) (r(dftrt_mideriv)) (`errcode')
		
		post `sim' (`r') ("MI compn main") (r(mcra)) (r(mderiv)) (r(beta0_mimain)) (r(betatrt_mimain)) ///
										       (r(se0_mimain)) (r(setrt_mimain)) (r(df0_mimain)) (r(dftrt_mimain)) (`errcode')
											   
		post `sim' (`r') ("MI compn by trt") (r(mcra)) (r(mderiv)) (r(beta0_mi2way)) (r(betatrt_mi2way)) ///
													   (r(se0_mi2way)) (r(setrt_mi2way)) (r(df0_mi2way)) (r(dftrt_mi2way)) (r(by1_err)) 
													   
		post `sim' (`r') ("MI compn by trt z1") (r(mcra)) (r(mderiv)) (r(beta0_mi3way)) (r(betatrt_mi3way)) ///
													      (r(se0_mi3way)) (r(setrt_mi3way)) (r(df0_mi3way)) (r(dftrt_mi3way)) (r(by2_err)) 			
	}
	 
	timer off 1	

}

timer list

postclose `sim'
postclose `postseed'

log close