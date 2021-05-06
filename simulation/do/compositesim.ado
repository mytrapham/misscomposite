//******************************************************************************
*Simulation for 'A comparison of methods for analysing a binary composite 
*endpoint with partially observed components in randomisedcontrolled trials'

*Tra My Pham | 04feb2020
//******************************************************************************
*Outcome types: 
*	- Simple: y = 1 if z1 = 1 or z2 = 1 or z3 = 1;
*	- Complex: y = 1 if z1 = 1 and (z2 = 1 or z3 = 1).

*Cases (related to lambda_123 in the log linear model for the components): 
*	- I (no 3-way interaction in either arm - noint); 
*	- II (3-way interaction in the control arm - intctr);
*	- III (3-way interaction in both arms - intrtctr).

*Missingness mechanisms of the components: 
*	- MCAR;
*	- MAR1 (conditional on randomised treatment x and the fully observed component z1);
* 	- MAR2 (conditional on randomised treatment x, the fully observed component z1, and their interaction).

*Methods for handling missing data: 
*	- CRA: 			complete record analysis;
*	- Deriv: 		analysis of the derived endpoint;
*	- MI-CRA: 		MI of the composite y after CRA;
*	- MI-Deriv: 	MI of composite y after Deriv; 
*	- MIC-main: 	MI of the components z2 and z3, x and z1 are included as main effects;
*	- MIC-x: 		MI of the components z2 and z3, z1 is included as main effect, imputation is stratified by x;
*	- MIC-x-z1: 	MI of the components z2 and z3, imputation is stratified by x and z1.

*Installation of -moremata- is required.
//******************************************************************************

program define compositesim, rclass
    version 15
	
	//Outcome types: 1-simple; 2-complex
	//Cases: 1-case I (noint); 2-case II (intctr); 3-case III (inttrtctr)
	//Missingness mechanisms: 1-MCAR; 2-MAR1; 3-MAR2
    syntax [ , outcometype(int 1) case(int 1) missingness(int 1) n(int 2000) m(int 5) burnin(int 5)] 
	
	qui	{
		set obs `n'
				
		//Generate randomised treatment
		gen x = runiform() <= 0.5
		
		//Define parameters for generating 8 combinations of the components 
		*Simple composite
		if `outcometype' == 1	{
			*Case I: noint
			if `case' == 1	{
				local lambdac = lambdac_casea_simp
				local lambdat = lambdat_casea_simp
				
				local lambdac_twoway = 2*`lambdac'+1
				local lambdac_threeway = 3*`lambdac'+3
				
				local lambdat_twoway = 2*`lambdat'+0.5
				local lambdat_threeway = 3*`lambdat'+1.5
			}
			else	{
				*Case II: intctr
				if `case' == 2	{
					local lambdac = lambdac_caseb_simp
					local lambdat = lambdat_caseb_simp
					
					local lambdac_twoway = 2*`lambdac'+1
					local lambdac_threeway = 3*`lambdac'+3.5
					
					local lambdat_twoway = 2*`lambdat'+0.5
					local lambdat_threeway = 3*`lambdat'+1.5
				}
				*Case III: inttrtctr
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
		*Complex composite
		else	{
				*Case I: noint
			if `case' == 1	{
				local lambdac = lambdac_casea_comp
				local lambdat = lambdat_casea_comp
				
				local lambdac_twoway = 2*`lambdac'+1
				local lambdac_threeway = 3*`lambdac'+3
				
				local lambdat_twoway = 2*`lambdat'+0.5
				local lambdat_threeway = 3*`lambdat'+1.5
			}
			else	{
				*Case II: intctr
				if `case' == 2	{
					local lambdac = lambdac_caseb_comp
					local lambdat = lambdat_caseb_comp
					
					local lambdac_twoway = 2*`lambdac'+1
					local lambdac_threeway = 3*`lambdac'+3.5
					
					local lambdat_twoway = 2*`lambdat'+0.5
					local lambdat_threeway = 3*`lambdat'+1.5		
				}
				*Case III: inttrtctr
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
		*Generate 8 combinations c of the components using the pre-selected values of lambdas
		
		local totc = exp(0) + 3*exp(`lambdac') + 3*exp(`lambdac_twoway') + exp(`lambdac_threeway')
		local tott = exp(0) + 3*exp(`lambdat') + 3*exp(`lambdat_twoway') + exp(`lambdat_threeway')
		
		gen c = .
		replace c = 1 if runiform() < exp(0)/`totc' & !x
		local totc = `totc' - exp(0)
		replace c = 2 if missing(c) & runiform() <= exp(`lambdac')/`totc' & !x
		local totc = `totc' - exp(`lambdac')
		replace c = 3 if missing(c) & runiform() <= exp(`lambdac')/`totc' & !x
		local totc = `totc' - exp(`lambdac')
		replace c = 4 if missing(c) & runiform() <= exp(`lambdac_twoway')/`totc' & !x 
		local totc = `totc' - exp(`lambdac_twoway')
		replace c = 5 if missing(c) & runiform() <= exp(`lambdac')/`totc' & !x
		local totc = `totc' - exp(`lambdac')
		replace c = 6 if missing(c) & runiform() <= exp(`lambdac_twoway')/`totc' & !x
		local totc = `totc' - exp(`lambdac_twoway')
		replace c = 7 if missing(c) & runiform() <= exp(`lambdac_twoway')/`totc' & !x
		replace c = 8 if missing(c) & !x 

		replace c = 1 if runiform() < exp(0)/`tott' & x
		local tott = `tott' - exp(0)
		replace c = 2 if missing(c) & runiform() <= exp(`lambdat')/`tott' & x
		local tott = `tott' - exp(`lambdat')
		replace c = 3 if missing(c) & runiform() <= exp(`lambdat')/`tott' & x
		local tott = `tott' - exp(`lambdat')
		replace c = 4 if missing(c) & runiform() <= exp(`lambdat_twoway')/`tott' & x 
		local tott = `tott' - exp(`lambdat_twoway')
		replace c = 5 if missing(c) & runiform() <= exp(`lambdat')/`tott' & x
		local tott = `tott' - exp(`lambdat')
		replace c = 6 if missing(c) & runiform() <= exp(`lambdat_twoway')/`tott' & x
		local tott = `tott' - exp(`lambdat_twoway')
		replace c = 7 if missing(c) & runiform() <= exp(`lambdat_twoway')/`tott' & x
		replace c = 8 if missing(c) & x 

		*Generate the components 
		gen z1 = c > 4
		gen z2 = inlist(c, 3, 4, 7, 8)
		gen z3 = !mod(c, 2) 

		//Generate the composite y
		*Simple composite
		if `outcometype' == 1	{
			gen y = inlist(1, z1, z2, z3)
		}
		*Complex composite
		else	{
			gen y = z1 & (z2 | z3)
		}
		
		//Specify parameter values for the selection models of the components
		*MCAR
		if `missingness' == 1	{
			local alpha0 = logit(0.7)	
			local alphax = 0
			local alphaz1 = 0
			local alphaxz1 = 0
		}
		else	{
				*MAR1
				if `missingness' == 2	{
					local alpha0 = 1.05
					local alphax = -0.75
					local alphaz1 = 0.25
					local alphaxz1 = 0
				}
				*MAR2
				else	{
					local alpha0 = 1.05
					local alphax = -0.75
					local alphaz1 = 0.25
					local alphaxz1 = 0.25
				}
		}
		
		//Introduce missing values in the components zs
		forval i = 2/3	{
			gen  z`i'mis = z`i' if runiform() <= invlogit(`alpha0' + `alphax'*x + `alphaz1'*z1 + `alphaxz1'*x*z1)
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
		logit y i.x
			return scalar beta0_fd = _b[_cons]
			return scalar betax_fd = _b[1.x]
			return scalar se0_fd = _se[_cons]
			return scalar sex_fd = _se[1.x]
			*For simsum 
			return scalar df0_fd = 1000
			return scalar dfx_fd = 1000
			
		
		//Complete record analysis
		logit ycra i.x
			return scalar beta0_cra = _b[_cons]
			return scalar betax_cra = _b[1.x]
			return scalar se0_cra = _se[_cons]
			return scalar sex_cra = _se[1.x]
			*For simsum 
			return scalar df0_cra = 1000
			return scalar dfx_cra = 1000

		//Derived endpoint
		logit yderiv i.x
			return scalar beta0_deriv = _b[_cons]
			return scalar betax_deriv = _b[1.x]
			return scalar se0_deriv = _se[_cons]
			return scalar sex_deriv = _se[1.x]
			*For simsum 
			return scalar df0_deriv = 1000
			return scalar dfx_deriv = 1000
			
		//MI at composite level (MI-CRA)	
		mi set wide
		mi register imputed ycra 
		mi register regular x 
		mi impute logit ycra i.x, add(`m') 
		mi estimate: logit ycra i.x
		
		mat def b = e(b_mi)
		mat def var = e(V_mi)
		mat def df = e(df_mi)

		return scalar beta0_micra  = b[1,3]
		return scalar betax_micra = b[1,2]	
		return scalar se0_micra = sqrt(var[3,3])
		return scalar sex_micra = sqrt(var[2,2])
		return scalar df0_micra = df[1,3]
		return scalar dfx_micra = df[1,2]
		
		//MI at composite level (MI-Deriv)
		mi extract 0, clear
		mi set wide
		mi register imputed yderiv 
		mi register regular x 
		mi impute logit yderiv i.x, add(`m') 

		mi estimate: logit yderiv i.x
		mat def b = e(b_mi)
		mat def var = e(V_mi)
		mat def df = e(df_mi)

		return scalar beta0_mideriv  = b[1,3]
		return scalar betax_mideriv = b[1,2]	
		return scalar se0_mideriv = sqrt(var[3,3])
		return scalar sex_mideriv = sqrt(var[2,2])
		return scalar df0_mideriv = df[1,3]
		return scalar dfx_mideriv = df[1,2]
		
		//MI at component level (treatment and z1 as main effects)
		mi extract 0, clear
		mi set wide
		mi register imputed z2mis z3mis
		mi register regular x z1
		mi impute chained (logit) z2mis z3mis = i.x i.z1, add(`m') burnin(`burnin')  
		
		*Simple composite
		if `outcometype' == 1	{
			mi passive: gen ycompon = inlist(1, z1, z2mis, z3mis)
		}
		*Complex composite
		else	{
			mi passive: gen ycompon = z1 & (z2mis | z3mis)
		}
		
		mi estimate: logit ycompon i.x
		mat def b = e(b_mi)
		mat def var = e(V_mi)
		mat def df = e(df_mi)

		return scalar beta0_micmain = b[1,3]
		return scalar betax_micmain = b[1,2]	
		return scalar se0_micmain = sqrt(var[3,3])
		return scalar sex_micmain = sqrt(var[2,2]) 
		return scalar df0_micmain = df[1,3]
		return scalar dfx_micmain = df[1,2]
		
		//MI at component level (z1 as main effect, stratified by treatment)
		mi extract 0, clear
		mi set wide
		mi register imputed z2mis z3mis
		mi register regular x z1
		
		cap mi impute chained (logit) z2mis z3mis = i.z1, by(x) add(`m') burnin(`burnin')  
		if _rc == 498	{
			return scalar by1_err = _rc
			mi impute chained (logit) z2mis z3mis = i.z1, by(x) add(`m') burnin(`burnin') augment 
		}
		else	{
			return scalar by1_err = .
		}
		
		*Simple composite
		if `outcometype' == 1	{
			mi passive: gen ycomponby = inlist(1, z1, z2mis, z3mis)
		}
		*Complex composite
		else	{
			mi passive: gen ycomponby = z1 & (z2mis | z3mis)
		}
		
		mi estimate: logit ycomponby i.x
		mat def b = e(b_mi)
		mat def var = e(V_mi)
		mat def df = e(df_mi)

		return scalar beta0_micx = b[1,3]
		return scalar betax_micx = b[1,2]	
		return scalar se0_micx = sqrt(var[3,3])
		return scalar sex_micx = sqrt(var[2,2])
		return scalar df0_micx = df[1,3]
		return scalar dfx_micx = df[1,2]

		//MI at component level (stratified by treatment and z1)
		cap drop ycomponby
		mi update	
		mi extract 0, clear
		mi set wide
		mi register imputed z2mis z3mis
		mi register regular x z1
		
		cap mi impute chained (logit) z2mis z3mis, by(z1 x) add(`m') burnin(`burnin')  
		if _rc == 498	{
			return scalar by2_err = _rc
			mi impute chained (logit) z2mis z3mis, by(z1 x) add(`m') burnin(`burnin') augment
		}
		else	{
			return scalar by2_err = .
		}
			
		*Simple composite	
		if `outcometype' == 1	{
			mi passive: gen ycomponby = inlist(1, z1, z2mis, z3mis)
		}
		*Complex composite
		else	{
			mi passive: gen ycomponby = z1 & (z2mis | z3mis)
		}
		
		mi estimate: logit ycomponby i.x
		mat def b = e(b_mi)
		mat def var = e(V_mi)
		mat def df = e(df_mi)

		return scalar beta0_micxz1 = b[1,3]
		return scalar betax_micxz1 = b[1,2]	
		return scalar se0_micxz1 = sqrt(var[3,3])
		return scalar sex_micxz1 = sqrt(var[2,2])
		return scalar df0_micxz1 = df[1,3]
		return scalar dfx_micxz1 = df[1,2]
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