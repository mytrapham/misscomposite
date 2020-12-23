//Outcome type: 1-simple; 2-complex
//Case: 1-case I (noint); 2-case II (intctr); 3-case III (inttrtctr)
//Missingness: 1-MCAR; 2-MAR1; 3-MAR2

*** SIMPLE INTCTR MAR2 ***
********************************************************************************

version 15

cd /lustre/scratch/scratch/zcakf79/Stata_output/

clear *
macro drop _all

cap log close
log using simple_intctr_mar2_2, smcl replace

set rng mt64 
set rngstream 2001
set seed 8022

local R = 1000

tempname sim
postfile `sim' int(rep) str28(method) float(mcra mderiv b0 btrt se0 setrt df0 dftrt) using simple_intctr_mar2_2, replace
	   			   
qui	{
	noi _dots 0, title("Simulation running...")
	timer on 1	

	forval r = 1/`R' 	{
		clear
		noi _dots `r' 0
		
		//Outcome type: 1-simple; 2-complex
		//Case: 1-case I (noint); 2-case II (intctr); 3-case III (inttrtctr)
		//Missingness: 1-MCAR; 2-MAR1; 3-MAR2
		compositesim, outcometype(1) case(2) missingness(3) n(2000) m(10) burnin(15)
								   
		post `sim' (`r') ("Full data") (r(mcra)) (r(mderiv)) (r(beta0_fd)) (r(betatrt_fd)) ///
									   (r(se0_fd)) (r(setrt_fd)) (r(df0_fd)) (r(dftrt_fd))
								 
		post `sim' (`r') ("CRA") (r(mcra)) (r(mderiv)) (r(beta0_cra)) (r(betatrt_cra)) ///
							     (r(se0_cra)) (r(setrt_cra)) (r(df0_cra)) (r(dftrt_cra))
		
		post `sim' (`r') ("Deriv") (r(mcra)) (r(mderiv)) (r(beta0_deriv)) (r(betatrt_deriv)) ///
								   (r(se0_deriv)) (r(setrt_deriv)) (r(df0_deriv)) (r(dftrt_deriv))
		
		post `sim' (`r') ("MI-CRA") (r(mcra)) (r(mderiv)) (r(beta0_micra)) (r(betatrt_micra)) ///
								    (r(se0_micra)) (r(setrt_micra)) (r(df0_micra)) (r(dftrt_micra))
		
		post `sim' (`r') ("MI-Deriv") (r(mcra)) (r(mderiv)) (r(beta0_mideriv)) (r(betatrt_mideriv)) ///
									  (r(se0_mideriv)) (r(setrt_mideriv)) (r(df0_mideriv)) (r(dftrt_mideriv))
		
		post `sim' (`r') ("MI compn main") (r(mcra)) (r(mderiv)) (r(beta0_mimain)) (r(betatrt_mimain)) ///
										       (r(se0_mimain)) (r(setrt_mimain)) (r(df0_mimain)) (r(dftrt_mimain))
											   
		post `sim' (`r') ("MI compn by trt") (r(mcra)) (r(mderiv)) (r(beta0_mi2way)) (r(betatrt_mi2way)) ///
													   (r(se0_mi2way)) (r(setrt_mi2way)) (r(df0_mi2way)) (r(dftrt_mi2way))	
													   
		post `sim' (`r') ("MI compn by trt z1") (r(mcra)) (r(mderiv)) (r(beta0_mi3way)) (r(betatrt_mi3way)) ///
													      (r(se0_mi3way)) (r(setrt_mi3way)) (r(df0_mi3way)) (r(dftrt_mi3way))		
	}
	 
	timer off 1	

}

timer list

postclose `sim'

log close