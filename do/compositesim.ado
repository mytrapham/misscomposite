//******************************************************************************
*Tra My Pham | 07dec2020
*Simulation for 'Analysis of a binary composite endpoint with partially observed 
*components in randomised controlled trials: a comparison of strategies'.

*Outcome type: 
*	- simple;
*	- complex.

*Case: 
*	- I (no 3-way interaction in either arm - noint); 
*	- II (3-way interaction in the control arm - intctr);
*	- III (3-way interaction in both arms - intrtrtctr).

*Missingness in the components: 
*	- MCAR;
*	- MAR1 (conditional on randomised treatment and the fully observed component);
* 	- MAR2 (conditional on randomised treatment, the fully observed component, and their interaction).

*Methods for missing data: 
*	- CRA;
*	- Deriv;
*	- MI-CRA;
*	- MI-Deriv; 
*	- MI compn main;
*	- MI compn by trt;
*	- MI compn by trt z1.
//******************************************************************************

program define compositesim, rclass
    version 15
	
	//Outcome type: 1-simple; 2-complex
	//Case: 1-case I (noint); 2-case II (intctr); 3-case III (inttrtctr)
	//Missingness: 1-MCAR; 2-MAR1; 3-MAR2
    syntax [ , outcometype(int 1) case(int 1) missingness(int 1) n(int 2000) m(int 5) burnin(int 5)] 
	
	qui	{
		set obs `n'
				
		//Generate randomised treatment
		gen trt = runiform() <= 0.5
		
		//Define parameters for generating 8 combinations of the components 
		if `outcometype' == 1	{
			*Case A: noint
			if `case' == 1	{
				local lambdac = lambdac_casea_simp
				local lambdat = lambdat_casea_simp
				
				local lambdac_twoway = 2*`lambdac'+1
				local lambdac_threeway = 3*`lambdac'+3
				
				local lambdat_twoway = 2*`lambdat'+0.5
				local lambdat_threeway = 3*`lambdat'+1.5
			}
			else	{
				*Case B: intctr
				if `case' == 2	{
					local lambdac = lambdac_caseb_simp
					local lambdat = lambdat_caseb_simp
					
					local lambdac_twoway = 2*`lambdac'+1
					local lambdac_threeway = 3*`lambdac'+3.5
					
					local lambdat_twoway = 2*`lambdat'+0.5
					local lambdat_threeway = 3*`lambdat'+1.5
				}
				*Case C: inttrtctr
				else	{
					local lambdac = lambdac_casec_simp
					local lambdat = lambdat_casec_simp
					
					local lambdac_twoway = `lambdac'
					local lambdac_threeway = `lambdac'
					
					local lambdat_twoway = `lambdat'
					local lambdat_threeway = `lambdat'					
				}
			}
		}
		
		else	{
				*Case A: noint
			if `case' == 1	{
				local lambdac = lambdac_casea_comp
				local lambdat = lambdat_casea_comp
				
				local lambdac_twoway = 2*`lambdac'+1
				local lambdac_threeway = 3*`lambdac'+3
				
				local lambdat_twoway = 2*`lambdat'+0.5
				local lambdat_threeway = 3*`lambdat'+1.5
			}
			else	{
				*Case B: intctr
				if `case' == 2	{
					local lambdac = lambdac_caseb_comp
					local lambdat = lambdat_caseb_comp
					
					local lambdac_twoway = 2*`lambdac'+1
					local lambdac_threeway = 3*`lambdac'+3.5
					
					local lambdat_twoway = 2*`lambdat'+0.5
					local lambdat_threeway = 3*`lambdat'+1.5		
				}
				*Case C: inttrtctr
				else	{
					local lambdac = lambdac_casec_comp
					local lambdat = lambdat_casec_comp
					
					local lambdac_twoway = 2*`lambdac'+1
					local lambdac_threeway = 3*`lambdac'+3.5
					
					local lambdat_twoway = 2*`lambdat'+0.5
					local lambdat_threeway = 3*`lambdat'+1
				}
			}
		} 

		//Generate the components zs
		*Generate 8 combinations c of the components using the selected lambdas
		
		local totc = exp(0) + 3*exp(`lambdac') + 3*exp(`lambdac_twoway') + exp(`lambdac_threeway')
		local tott = exp(0) + 3*exp(`lambdat') + 3*exp(`lambdat_twoway') + exp(`lambdat_threeway')
		
		gen c = .
		replace c = 1 if runiform() < exp(0)/`totc' & !trt
		local totc = `totc' - exp(0)
		replace c = 2 if missing(c) & runiform() <= exp(`lambdac')/`totc' & !trt
		local totc = `totc' - exp(`lambdac')
		replace c = 3 if missing(c) & runiform() <= exp(`lambdac')/`totc' & !trt
		local totc = `totc' - exp(`lambdac')
		replace c = 4 if missing(c) & runiform() <= exp(`lambdac_twoway')/`totc' & !trt 
		local totc = `totc' - exp(`lambdac_twoway')
		replace c = 5 if missing(c) & runiform() <= exp(`lambdac')/`totc' & !trt
		local totc = `totc' - exp(`lambdac')
		replace c = 6 if missing(c) & runiform() <= exp(`lambdac_twoway')/`totc' & !trt
		local totc = `totc' - exp(`lambdac_twoway')
		replace c = 7 if missing(c) & runiform() <= exp(`lambdac_twoway')/`totc' & !trt
		replace c = 8 if missing(c) & !trt 

		replace c = 1 if runiform() < exp(0)/`tott' & trt
		local tott = `tott' - exp(0)
		replace c = 2 if missing(c) & runiform() <= exp(`lambdat')/`tott' & trt
		local tott = `tott' - exp(`lambdat')
		replace c = 3 if missing(c) & runiform() <= exp(`lambdat')/`tott' & trt
		local tott = `tott' - exp(`lambdat')
		replace c = 4 if missing(c) & runiform() <= exp(`lambdat_twoway')/`tott' & trt 
		local tott = `tott' - exp(`lambdat_twoway')
		replace c = 5 if missing(c) & runiform() <= exp(`lambdat')/`tott' & trt
		local tott = `tott' - exp(`lambdat')
		replace c = 6 if missing(c) & runiform() <= exp(`lambdat_twoway')/`tott' & trt
		local tott = `tott' - exp(`lambdat_twoway')
		replace c = 7 if missing(c) & runiform() <= exp(`lambdat_twoway')/`tott' & trt
		replace c = 8 if missing(c) & trt 

		*Generate the components 
		gen z1 = c > 4
		gen z2 = inlist(c, 3, 4, 7, 8)
		gen z3 = !mod(c, 2) 

		//Generate the composite y
		if `outcometype' == 1	{
			gen y = inlist(1, z1, z2, z3)
		}
		else	{
			gen y = z1 & (z2 | z3)
		}
		
		//Generate partially observed components
		*MCAR
		if `missingness' == 1	{
			local alpha0 = logit(0.7)	
			local alphax = 0
			local alphaz1 = 0
			local alphaxz = 0
		}
		else	{
				*MAR 1
				if `missingness' == 2	{
					local alpha0 = 1.05
					local alphax = -0.75
					local alphaz1 = 0.25
					local alphaxz = 0
				}
				*MAR 2
				else	{
					local alpha0 = 1.05
					local alphax = -0.75
					local alphaz1 = 0.25
					local alphaxz = 0.25
				}
		}
		
		//Introduce missing values in the components zs
		forval i = 2/3	{
			gen  z`i'mis = z`i' if runiform() <= invlogit(`alpha0' + `alphax'*trt + `alphaz1'*z1 + `alphaxz'*trt*z1)
			mark r`i' if !missing(z`i'mis)
		}
		
		//Mark complete records
		gen ycra = y if !missing(z2mis) & !missing(z3mis)
		count if missing(ycra)
		return scalar mcra = r(N)/_N
		
		//Mark derived endpoint
		gen yderiv = ycra
		if `outcometype' == 1	{
			replace yderiv = 1 if inlist(1, z1, z2mis, z3mis) & missing(yderiv)
		}
		else	{
			replace yderiv = 1 if z1 & (z2mis | z3mis) & missing(yderiv)
		}
		count if missing(yderiv)
		return scalar mderiv = r(N)/_N
		
		//Full data 
		logit y i.trt
			return scalar beta0_fd = _b[_cons]
			return scalar betatrt_fd = _b[1.trt]
			return scalar se0_fd = _se[_cons]
			return scalar setrt_fd = _se[1.trt]
			*For simsum 
			return scalar df0_fd = 1000
			return scalar dftrt_fd = 1000
			
		
		//Complete record analysis
		logit ycra i.trt
			return scalar beta0_cra = _b[_cons]
			return scalar betatrt_cra = _b[1.trt]
			return scalar se0_cra = _se[_cons]
			return scalar setrt_cra = _se[1.trt]
			*For simsum 
			return scalar df0_cra = 1000
			return scalar dftrt_cra = 1000

		//Derived endpoint
		logit yderiv i.trt
			return scalar beta0_deriv = _b[_cons]
			return scalar betatrt_deriv = _b[1.trt]
			return scalar se0_deriv = _se[_cons]
			return scalar setrt_deriv = _se[1.trt]
			*For simsum 
			return scalar df0_deriv = 1000
			return scalar dftrt_deriv = 1000
			
		//MI at composite level (MI-CRA)	
		mi set wide
		mi register imputed ycra 
		mi register regular trt 
		mi impute logit ycra i.trt, add(`m') 
		mi estimate: logit ycra i.trt
		
		mat def b = e(b_mi)
		mat def var = e(V_mi)
		mat def df = e(df_mi)

		return scalar beta0_micra  = b[1,3]
		return scalar betatrt_micra = b[1,2]	
		return scalar se0_micra = sqrt(var[3,3])
		return scalar setrt_micra = sqrt(var[2,2])
		return scalar df0_micra = df[1,3]
		return scalar dftrt_micra = df[1,2]
		
		//MI at composite level (MI-Deriv)
		mi extract 0, clear
		mi set wide
		mi register imputed yderiv 
		mi register regular trt 
		mi impute logit yderiv i.trt, add(`m') 

		mi estimate: logit yderiv i.trt
		mat def b = e(b_mi)
		mat def var = e(V_mi)
		mat def df = e(df_mi)

		return scalar beta0_mideriv  = b[1,3]
		return scalar betatrt_mideriv = b[1,2]	
		return scalar se0_mideriv = sqrt(var[3,3])
		return scalar setrt_mideriv = sqrt(var[2,2])
		return scalar df0_mideriv = df[1,3]
		return scalar dftrt_mideriv = df[1,2]
		
		//MI at component level (treatment as main effect)
		mi extract 0, clear
		mi set wide
		mi register imputed z2mis z3mis
		mi register regular trt z1
		mi impute chained (logit) z2mis z3mis = i.trt i.z1, add(`m') burnin(`burnin')  
		if `outcometype' == 1	{
			mi passive: gen ycompon = inlist(1, z1, z2mis, z3mis)
		}
		else	{
			mi passive: gen ycompon = z1 & (z2mis | z3mis)

		}
		mi estimate: logit ycompon i.trt
		mat def b = e(b_mi)
		mat def var = e(V_mi)
		mat def df = e(df_mi)

		return scalar beta0_mimain = b[1,3]
		return scalar betatrt_mimain = b[1,2]	
		return scalar se0_mimain = sqrt(var[3,3])
		return scalar setrt_mimain = sqrt(var[2,2])
		return scalar df0_mimain = df[1,3]
		return scalar dftrt_mimain = df[1,2]
		
		//MI at component level (by treatment)
		mi extract 0, clear
		mi set wide
		mi register imputed z2mis z3mis
		mi register regular trt z1
		
		cap mi impute chained (logit) z2mis z3mis = i.z1, by(trt) add(`m') burnin(`burnin')  
		if _rc == 498	{
			return scalar by1_err = _rc
			mi impute chained (logit) z2mis z3mis = i.z1, by(trt) add(`m') burnin(`burnin') augment 
		}
		else	{
			return scalar by1_err = .
		}
			
		if `outcometype' == 1	{
			mi passive: gen ycomponby = inlist(1, z1, z2mis, z3mis)
		}
		else	{
			mi passive: gen ycomponby = z1 & (z2mis | z3mis)
		}
		
		mi estimate: logit ycomponby i.trt
		mat def b = e(b_mi)
		mat def var = e(V_mi)
		mat def df = e(df_mi)

		return scalar beta0_mi2way = b[1,3]
		return scalar betatrt_mi2way = b[1,2]	
		return scalar se0_mi2way = sqrt(var[3,3])
		return scalar setrt_mi2way = sqrt(var[2,2])
		return scalar df0_mi2way = df[1,3]
		return scalar dftrt_mi2way = df[1,2]

		//MI at component level (by treatment and z1)
		cap drop ycomponby
		mi update	
		mi extract 0, clear
		mi set wide
		mi register imputed z2mis z3mis
		mi register regular trt z1
		
		cap mi impute chained (logit) z2mis z3mis, by(z1 trt) add(`m') burnin(`burnin')  
		if _rc == 498	{
			return scalar by2_err = _rc
			mi impute chained (logit) z2mis z3mis, by(z1 trt) add(`m') burnin(`burnin') augment
		}
		else	{
			return scalar by2_err = .
		}
			
		if `outcometype' == 1	{
			mi passive: gen ycomponby = inlist(1, z1, z2mis, z3mis)
		}
		else	{
			mi passive: gen ycomponby = z1 & (z2mis | z3mis)

		}
		mi estimate: logit ycomponby i.trt
		mat def b = e(b_mi)
		mat def var = e(V_mi)
		mat def df = e(df_mi)

		return scalar beta0_mi3way = b[1,3]
		return scalar betatrt_mi3way = b[1,2]	
		return scalar se0_mi3way = sqrt(var[3,3])
		return scalar setrt_mi3way = sqrt(var[2,2])
		return scalar df0_mi3way = df[1,3]
		return scalar dftrt_mi3way = df[1,2]
	}
	
end





//Calculate paramters of the log linear models used for generating the components
mata

	//// SIMPLE COMPOSITE
	a_c = ((exp(0.3)/(1+exp(0.3)))) * (1/(1-(exp(0.3)/(1+exp(0.3)))))*1/3
	a_t = ((exp(1.65)/(1+exp(1.65)))) * (1/(1-(exp(1.65)/(1+exp(1.65)))))*1/3
	
	////////////////////////////////////////////////////////////////////////////
	//Case I: noint 
	////////////////////////////////////////////////////////////////////////////
	/*Assume lambda_1 = lambda_2 = lambda_3 = lambdac; lambda_12 = lambda_13 = lambda_23 = 1; lambda_123 = 0
	*Calculate lambdac for the control arm*/ 
	function myfuncc_casea_simp(x_c, a_c) return(exp(x_c) + exp(2*x_c+1) + (1/3)*exp(3*x_c+3) - a_c)
	mm_root(x_c=., &myfuncc_casea_simp(), -3, 3, 0.0001, 100000, a_c)
	st_numscalar("lambdac_casea_simp", x_c)

	/*Assume lambda_1 = lambda_2 = lambda_3 = lambdat; lambda_12 = lambda_13 = lambda_23 = 0.5; lambda_123 = 0
	*Calculate lambdat for the treatment arm*/ 
	function myfunct_casea_simp(x_t, a_t) return(exp(x_t) + exp(2*x_t+0.5) + (1/3)*exp(3*x_t+1.5) - a_t)
	mm_root(x_t=., &myfunct_casea_simp(), -3, 3, 0.0001, 100000, a_t)
	st_numscalar("lambdat_casea_simp", x_t)
	
	
	////////////////////////////////////////////////////////////////////////////
	//Case II: intctr 
	////////////////////////////////////////////////////////////////////////////
	/*Assume lambda_1 = lambda_2 = lambda_3 = lambdac; lambda_12 = lambda_13 = lambda_23 = 1; lambda_123 = 0.5
	*Calculate lambdac for the control arm*/ 
	function myfuncc_caseb_simp(x_c, a_c) return(exp(x_c) + exp(2*x_c+1) + (1/3)*exp(3*x_c+3.5) - a_c)
	mm_root(x_c=., &myfuncc_caseb_simp(), -3, 3, 0.0001, 100000, a_c)
	st_numscalar("lambdac_caseb_simp", x_c)

	/*Assume lambda_1 = lambda_2 = lambda_3 = lambdat; lambda_12 = lambda_13 = lambda_23 = 0.5; lambda_123 = 0
	*Calculate lambdat for the treatment arm*/ 
	function myfunct_caseb_simp(x_t, a_t) return(exp(x_t) + exp(2*x_t+0.5) + (1/3)*exp(3*x_t+1.5) - a_t)
	mm_root(x_t=., &myfunct_caseb_simp(), -3, 3, 0.0001, 100000, a_t)
	st_numscalar("lambdat_caseb_simp", x_t)

	
	////////////////////////////////////////////////////////////////////////////
	//Case III: inttrtctr 
	////////////////////////////////////////////////////////////////////////////
	/*Assume lambda_1 = lambda_2 = lambda_3 = lambdac; lambda_12 = lambda_13 = lambda_23 = -lambdac; lambda_123 = lambdac
	*Calculate lambdat for the treatment arm*/ 
	x_c = ln(a_c/(7/3))
	st_numscalar("lambdac_casec_simp", x_c)

	/*Assume lambda_1 = lambda_2 = lambda_3 = lambdat; lambda_12 = lambda_13 = lambda_23 = -lambdat; lambda_123 = lambdat
	*Calculate lambdat for the treatment arm*/ 	
	x_t = ln(a_t/(7/3))
	st_numscalar("lambdat_casec_simp", x_t)

	
	
	
	
	//// COMPLEX COMPOSITE
	b_c = 1 - (exp(0.3)/(1+exp(0.3)))
	b_t = 1 - (exp(1.65)/(1+exp(1.65)))

	////////////////////////////////////////////////////////////////////////////
	//Case I: noint 
	////////////////////////////////////////////////////////////////////////////
	/*Assume lambda_1 = lambda_2 = lambda_3 = lambdac; lambda_12 = lambda_13 = lambda_23 = 1; lambda_123 = 0
	*Calculate lambdac for the control arm*/ 
	function myfuncc_casea_comp(x_c, b_c) return((3-3*b_c)*exp(x_c) + (1-3*b_c)*exp(2*x_c+1) - b_c*exp(3*x_c+3) + (1-b_c))
	mm_root(x_c=., &myfuncc_casea_comp(), -3, 3, 0.0001, 100000, b_c)
	st_numscalar("lambdac_casea_comp", x_c)
	
	/*Assume lambda_1 = lambda_2 = lambda_3 = lambdat; lambda_12 = lambda_13 = lambda_23 = 0.5; lambda_123 = 0
	*Calculate lambdat for the treatment arm*/ 
	function myfunct_casea_comp(x_t, b_t) return((3-3*b_t)*exp(x_t) + (1-3*b_t)*exp(2*x_t+0.5) - b_t*exp(3*x_t+1.5) + (1-b_t))
	mm_root(x_t=., &myfunct_casea_comp(), -3, 3, 0.0001, 100000, b_t)
	st_numscalar("lambdat_casea_comp", x_t)
	
	
	////////////////////////////////////////////////////////////////////////////
	//Case II: intctr 
	////////////////////////////////////////////////////////////////////////////
	/*Assume lambda_1 = lambda_2 = lambda_3 = lambdac; lambda_12 = lambda_13 = lambda_23 = 1; lambda_123 = 0.5
	*Calculate lambdac for the control arm*/ 
	function myfuncc_caseb_comp(x_c, b_c) return((3-3*b_c)*exp(x_c) + (1-3*b_c)*exp(2*x_c+1) - b_c*exp(3*x_c+3.5) + (1-b_c))
	mm_root(x_c=., &myfuncc_caseb_comp(), -3, 3, 0.0001, 100000, b_c)
	st_numscalar("lambdac_caseb_comp", x_c)

	/*Assume lambda_1 = lambda_2 = lambda_3 = lambdat; lambda_12 = lambda_13 = lambda_23 = 0.5; lambda_123 = 0
	*Calculate lambdat for the treatment arm*/ 
	function myfunct_caseb_comp(x_t, b_t) return((3-3*b_t)*exp(x_t) + (1-3*b_t)*exp(2*x_t+0.5) - b_t*exp(3*x_t+1.5) + (1-b_t))
	mm_root(x_t=., &myfunct_caseb_comp(), -3, 3, 0.0001, 100000, b_t)
	st_numscalar("lambdat_caseb_comp", x_t)
	
	
	////////////////////////////////////////////////////////////////////////////
	//Case III: inttrtctr 
	////////////////////////////////////////////////////////////////////////////
	/*Assume lambda_1 = lambda_2 = lambda_3 = lambdac; lambda_12 = lambda_13 = lambda_23 = 1; lambda_123 = 0.5
	*Calculate lambdac for the control arm*/ 
	function myfuncc_casec_comp(x_c, b_c) return((3-3*b_c)*exp(x_c) + (1-3*b_c)*exp(2*x_c+1) - b_c*exp(3*x_c+3.5) + (1-b_c))
	mm_root(x_c=., &myfuncc_casec_comp(), -3, 3, 0.0001, 100000, b_c)
	st_numscalar("lambdac_casec_comp", x_c)

	/*Assume lambda_1 = lambda_2 = lambda_3 = lambdat; lambda_12 = lambda_13 = lambda_23 = 0.5; lambda_123 = -0.5
	*Calculate lambdat for the treatment arm*/ 
	function myfunct_casec_comp(x_t, b_t) return((3-3*b_t)*exp(x_t) + (1-3*b_t)*exp(2*x_t+0.5) - b_t*exp(3*x_t+1) + (1-b_t))
	mm_root(x_t=., &myfunct_casec_comp(), -3, 3, 0.0001, 100000, b_t)
	st_numscalar("lambdat_casec_comp", x_t)
	
end