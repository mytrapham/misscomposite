////////////////////////////////////////////////////////////////////////////////
*Re-analysis of TOPPS
*Authors: T M Pham & B C Kahan
*Date: 15dec2020
********************************************************************************
*Notes:
/*
Variables in original TOPPS data
bleedgrade1-30 (bleed grade 0-4, days 1-30 (grade 2-4 = outcome event))
whograde2 (whether patient had 1+ grade 2-4 bleeds during follow-up)
treat (treatment allocation)
diagnosis (minimisation var)
treatment_plan (minimisation var)
trialno (unique trial id)

Non-convergence occured for MI compn by trt implemented by -mi impute chained-; 
alternatives considered include:
	- Using -ice- instead of -mi impute chained-;
	- Conditioning on 2 neighbouring blocks in -mi impute chained-
*/
////////////////////////////////////////////////////////////////////////////////


version 15

ssc install ice

clear *

set rng mt64 
set seed 21758

//Define paths
local datapath "C:\Users\wew584\Documents\Work\2_Studies\TOPPS\Data"
local outputpath "C:\Users\wew584\Dropbox\Methodology\9_16 Missing composites (Tra)"

//Log output
cap log close
log using "`outputpath'\topps_reanalysis_log", smcl replace

//Save random number states for re-creating MI results
tempname postseed
postfile `postseed' int(approach) str28(method) str2000(s1 s2) str1100(s3) using "`outputpath'\topps_reanalysis_miseed.dta", replace

//Save re-analysis results
tempname topps
postfile `topps' int(approach) str28(method) str9(param) float(b_trt ll_trt ul_trt p) ///
				 using "`outputpath'\topps_reanalysis.dta", replace
	
//Load original TOPPS data
*For the purpose of this re-analysis we will use bleedgrade1-30, treat, trialno
cd "`outputpath'"
use "`datapath'\TOPPS - cleaned.dta", clear

// Generate indicator of whether each participant had a grade 2-4 bleed each day during follow-up (0, 1, .)
rename whograde2 orig_whograde2
rename whograde3 orig_whograde3

forval i = 1/30	{
	gen whograde`i' = cond(missing(bleedgrade`i'), ., cond(bleedgrade`i' < 2, 0, 1))

}

misstable patterns whograde1-whograde30, freq

keep trialno treat bleedgrade* whograde*
order trialno treat bleedgrade* whograde*

compress

save "`outputpath'\topps_data_formi.dta", replace






	/*Approach 1: 
	- Each block is set to missing if bleeding status is missing for at least 1 day
	- Each block takes value 1 if bleeding occurs on at least 1 day
	*/

use "`outputpath'\topps_data_formi.dta", clear

//Generate 6 blocks of 5 days
local i = 1
local j = 5
forval b = 1/6	{
	tempvar block`b'_miss
	egen block`b' = rowmax(whograde`i'-whograde`j')
		
	*Check number of days with missing bleeding status within each block (overall and by treatment)
	egen `block`b'_miss' = rowmiss(whograde`i'-whograde`j')
	tab `block`b'_miss' treat, col
	
	replace block`b' = . if `block`b'_miss' > 0
		
	*Check number of participants with at least 1 bleeding episode in each block (overall and by treatment)
	tab block`b' treat, m col
	
	local i = `i' + 5
	local j = `j' + 5
}

misstable patterns block1-block6, freq

//Complete record analysis - composite y with 6 components
tempvar cra_miss
egen `cra_miss' = rowmiss(block1-block6)
gen ycra = 1 if inlist(1, block1, block2, block3, block4, block5, block6) & `cra_miss' == 0
replace ycra = 0 if missing(ycra) & `cra_miss' == 0

*Check % missing values in the composite and % of events by treatment
tempvar r_ycra
mark `r_ycra' if !missing(ycra)
tab `r_ycra' treat, col
tab ycra treat, col

*Substantive analysis model - risk difference
glm ycra i.treat, link(identity) family(binomial) iterate(50)
	mat def R 		= r(table)
	local b_trt 	= R[1, 2]
	local ll_trt 	= R[5, 2]
	local ul_trt 	= R[6, 2]
	local p 		= R[4, 2]

post `topps' (1) ("CRA") ("Risk diff") (`b_trt') (`ll_trt') (`ul_trt') (`p')

*Substantive analysis model - log OR
logit ycra i.treat, iterate(50)
	mat def R 		= r(table)
	local b_trt 	= R[1, 2]
	local ll_trt 	= R[5, 2]
	local ul_trt 	= R[6, 2]
	local p 		= R[4, 2]

post `topps' (1) ("CRA") ("Log OR") (`b_trt') (`ll_trt') (`ul_trt') (`p')

//Derived composite
gen yderiv = ycra
replace yderiv = 1 if inlist(1, block1, block2, block3, block4, block5, block6) & missing(yderiv)

*Check % missing values in the composite and % of events by treatment
tempvar r_yderiv
mark `r_yderiv' if !missing(yderiv)
tab `r_yderiv' treat, col
tab yderiv treat, col

*Substative analysis model - risk difference
glm yderiv i.treat, link(identity) family(binomial) iterate(50)
	mat def R 		= r(table)
	local b_trt 	= R[1, 2]
	local ll_trt 	= R[5, 2]
	local ul_trt 	= R[6, 2]
	local p 		= R[4, 2]

post `topps' (1) ("Deriv") ("Risk diff") (`b_trt') (`ll_trt') (`ul_trt') (`p')

*Substantive analysis model - log OR
logit yderiv i.treat, iterate(50)
	mat def R 		= r(table)
	local b_trt 	= R[1, 2]
	local ll_trt 	= R[5, 2]
	local ul_trt 	= R[6, 2]
	local p 		= R[4, 2]

post `topps' (1) ("Deriv") ("Log OR") (`b_trt') (`ll_trt') (`ul_trt') (`p')

//MI at composite level (MI-CRA)
post `postseed' (1) ("MI-CRA") (substr(c(rngstate),1,2000)) (substr(c(rngstate),2001,2000)) (substr(c(rngstate),4001,.))

local m = 50
local c = 20
compress

mi set wide
mi register imputed ycra 
mi register regular treat 

mi impute logit ycra i.treat, add(`m') augment iterate(50) dots

*Check imputed values in the first 5 imputed datasets
forval i = 1/5	{
	tab _`i'_ycra _mi_m, col
}

*Check % of events in the composite by treatment 
mi estimate: prop ycra, over(treat)

*Substantive analysis model - risk difference
mi estimate, mcerr vart: glm ycra i.treat, link(identity) family(binomial) iterate(50)
	mat def R 		= r(table)
	local b_trt 	= R[1, 2]
	local ll_trt 	= R[5, 2]
	local ul_trt 	= R[6, 2]
	local p 		= R[4, 2]

post `topps' (1) ("MI-CRA") ("Risk diff") (`b_trt') (`ll_trt') (`ul_trt') (`p')

*Substantive analysis model - log OR
mi estimate, mcerr vart: logit ycra i.treat, iterate(50)
	mat def R 		= r(table)
	local b_trt 	= R[1, 2]
	local ll_trt 	= R[5, 2]
	local ul_trt 	= R[6, 2]
	local p 		= R[4, 2]

post `topps' (1) ("MI-CRA") ("Log OR") (`b_trt') (`ll_trt') (`ul_trt') (`p')

//MI at composite level (MI-Deriv)
post `postseed' (1) ("MI-Deriv") (substr(c(rngstate),1,2000)) (substr(c(rngstate),2001,2000)) (substr(c(rngstate),4001,.))

mi extract 0, clear
mi set wide
mi register imputed yderiv 
mi register regular treat 

mi impute logit yderiv i.treat, add(`m') augment iterate(50) dots

*Check imputed values in the first 5 imputed datasets
forval i = 1/5	{
	tab _`i'_yderiv _mi_m, col
}

*Check % of events in the composite by treatment 
mi estimate: prop yderiv, over(treat)

*Substantive analysis model - risk difference
mi estimate, mcerr vart: glm yderiv i.treat, link(identity) family(binomial) iterate(50)
	mat def R 		= r(table)
	local b_trt 	= R[1, 2]
	local ll_trt 	= R[5, 2]
	local ul_trt 	= R[6, 2]
	local p 		= R[4, 2]

post `topps' (1) ("MI-Deriv") ("Risk diff") (`b_trt') (`ll_trt') (`ul_trt') (`p')

*Substantive analysis model - log OR
mi estimate, mcerr vart: logit yderiv i.treat, iterate(50)
	mat def R 		= r(table)
	local b_trt 	= R[1, 2]
	local ll_trt 	= R[5, 2]
	local ul_trt 	= R[6, 2]
	local p 		= R[4, 2]

post `topps' (1) ("MI-Deriv") ("Log OR") (`b_trt') (`ll_trt') (`ul_trt') (`p')

//MI at component level (MI compn main)
post `postseed' (1) ("MI compn main") (substr(c(rngstate),1,2000)) (substr(c(rngstate),2001,2000)) (substr(c(rngstate),4001,.))

mi extract 0, clear
mi set wide
mi register imputed block1-block6
mi register regular treat 

mi impute chained (logit, augment iterate(50)) block1 block2 block3 block4 block5 block6 = i.treat, add(`m') burnin(`c') dots 

*Check imputed values in the e.g. 5th imputed dataset
forval i = 1/6	{
	tab _5_block`i' _mi_m, col
}

*Check % of events in each block by treatment
mi estimate: prop block1-block6, over(treat)

*Generate composite 
mi passive: gen ycompnmain = inlist(1, block1, block2, block3, block4, block5, block6)

*Check passively imputed values in the first 5 imputed datasets
forval i = 1/5	{
	tab _`i'_ycompnmain _mi_m, col
}

*Check % of events in the composite by treatment 
mi estimate: prop ycompnmain, over(treat)

*Substantive analysis model - risk difference
mi estimate, mcerr vart: glm ycompnmain i.treat, link(identity) family(binomial) iterate(50)
	mat def R 		= r(table)
	local b_trt 	= R[1, 2]
	local ll_trt 	= R[5, 2]
	local ul_trt 	= R[6, 2]
	local p 		= R[4, 2]

post `topps' (1) ("MI compn main") ("Risk diff") (`b_trt') (`ll_trt') (`ul_trt') (`p')

*Substantive analysis model - log OR
mi estimate, mcerr vart: logit ycompnmain i.treat, iterate(50)
	mat def R 		= r(table)
	local b_trt 	= R[1, 2]
	local ll_trt 	= R[5, 2]
	local ul_trt 	= R[6, 2]
	local p 		= R[4, 2]

post `topps' (1) ("MI compn main") ("Log OR") (`b_trt') (`ll_trt') (`ul_trt') (`p')

cap drop ycompnmain
mi update

//MI at component level (MI compn by trt)

*MI using -mi impute chained-, conditioning on all other blocks when imputing each block (failed due to non-convergence in past runs)
post `postseed' (1) ("MI compn by trt") (substr(c(rngstate),1,2000)) (substr(c(rngstate),2001,2000)) (substr(c(rngstate),4001,.))

mi extract 0, clear
mi set wide
mi register imputed block1-block6
mi register regular treat 

cap mi impute chained (logit, augment iterate(50)) block1 block2 block3 block4 block5 block6, by(treat) add(`m') burnin(`c')  

if !_rc	{		  
	*Check imputed values in the e.g. 5th imputed dataset
	forval i = 1/6	{
		tab _5_block`i' _mi_m, col
	}

	*Check % of events in each block by treatment
	mi estimate: prop block1-block6, over(treat)

	*Generate composite
	mi passive: gen ycompnby = inlist(1, block1, block2, block3, block4, block5, block6)

	*Check passively imputed values in the first 5 imputed datasets
	forval i = 1/5	{
		tab _`i'_ycompnby _mi_m, col
	}

	*Check % of events in the composite by treatment 
	mi estimate: prop ycompnby, over(treat)

	*Substantive analysis model - risk difference
	mi estimate, mcerr vart: glm ycompnby i.treat, link(identity) family(binomial) iterate(50)
		mat def R 		= r(table)
		local b_trt 	= R[1, 2]
		local ll_trt 	= R[5, 2]
		local ul_trt 	= R[6, 2]
		local p 		= R[4, 2]

	post `topps' (1) ("MI compn by trt") ("Risk diff") (`b_trt') (`ll_trt') (`ul_trt') (`p')

	*Substantive analysis model - log OR
	mi estimate, mcerr vart: logit ycompnby i.treat, iterate(50)
		mat def R 		= r(table)
		local b_trt 	= R[1, 2]
		local ll_trt 	= R[5, 2]
		local ul_trt 	= R[6, 2]
		local p 		= R[4, 2]

	post `topps' (1) ("MI compn by trt") ("Log OR") (`b_trt') (`ll_trt') (`ul_trt') (`p')
	
	cap drop ycompnby
	mi update	
}

else	{
	di as err "MI compn by trt failed with error code " as result _rc
}
	
*MI using -mi impute chained-, conditioning on 2 neighbouring blocks when imputing each block
post `postseed' (1) ("MI compn by trt cond2") (substr(c(rngstate),1,2000)) (substr(c(rngstate),2001,2000)) (substr(c(rngstate),4001,.))

mi extract 0, clear
mi set wide
mi register imputed block1-block6
mi register regular treat

cap mi impute chained (logit, omit(i.block4 i.block5 i.block6) augment iterate(50)) block1 ///
					  (logit, omit(i.block4 i.block5 i.block6) augment iterate(50)) block2 ///
				      (logit, omit(i.block1 i.block5 i.block6) augment iterate(50)) block3 ///
					  (logit, omit(i.block1 i.block2 i.block6) augment iterate(50)) block4 ///
				      (logit, omit(i.block1 i.block2 i.block3) augment iterate(50)) block5 ///
					  (logit, omit(i.block1 i.block2 i.block3) augment iterate(50)) block6 ///
				      , by(treat) add(`m') burnin(`c')
					  
if !_rc	{		  
	*Check imputed values in the e.g. 5th imputed dataset
	forval i = 1/6	{
		tab _5_block`i' _mi_m, col
	}

	*Check % of events in each block by treatment
	mi estimate: prop block1-block6, over(treat)

	*Generate composite
	mi passive: gen ycompnby = inlist(1, block1, block2, block3, block4, block5, block6)

	*Check passively imputed values in the first 5 imputed datasets
	forval i = 1/5	{
		tab _`i'_ycompnby _mi_m, col
	}

	*Check % of events in the composite by treatment 
	mi estimate: prop ycompnby, over(treat)

	*Substantive analysis model - risk difference
	mi estimate, mcerr vart: glm ycompnby i.treat, link(identity) family(binomial) iterate(50)
		mat def R 		= r(table)
		local b_trt 	= R[1, 2]
		local ll_trt 	= R[5, 2]
		local ul_trt 	= R[6, 2]
		local p 		= R[4, 2]

	post `topps' (1) ("MI compn by trt cond2") ("Risk diff") (`b_trt') (`ll_trt') (`ul_trt') (`p')

	*Substantive analysis model - log OR
	mi estimate, mcerr vart: logit ycompnby i.treat, iterate(50)
		mat def R 		= r(table)
		local b_trt 	= R[1, 2]
		local ll_trt 	= R[5, 2]
		local ul_trt 	= R[6, 2]
		local p 		= R[4, 2]

	post `topps' (1) ("MI compn by trt cond2") ("Log OR") (`b_trt') (`ll_trt') (`ul_trt') (`p')
	
	cap drop ycompnby
	mi update
}

else	{
	di as err "MI compn by trt cond2 failed with error code " as result _rc
}

*MI using -ice-, conditioning on all other blocks when imputing each block (uvis automatically checks for perfect prediction)
post `postseed' (1) ("MI compn by trt ice") (substr(c(rngstate),1,2000)) (substr(c(rngstate),2001,2000)) (substr(c(rngstate),4001,.))

mi extract 0, clear

cap ice block1 block2 block3 block4 block5 block6, by(treat) m(`m') cycles(`c') saving("`outputpath'\mi_compn_by_trt_ice_approach1.dta", replace) 

if !_rc	{
	use "`outputpath'\mi_compn_by_trt_ice_approach1.dta", clear

	mi import ice, imputed(block1-block6) clear
	mi convert wide, clear
													 
	*Check imputed values in the e.g. 5th imputed dataset
	forval i = 1/6	{
		tab _5_block`i' _mi_m, col
	}

	*Check % of events in each block by treatment
	mi estimate: prop block1-block6, over(treat)

	*Generate composite
	mi passive: gen ycompnby = inlist(1, block1, block2, block3, block4, block5, block6)

	*Check passively imputed values in the first 5 imputed datasets
	forval i = 1/5	{
		tab _`i'_ycompnby _mi_m, col
	}

	*Check % of events in the composite by treatment 
	mi estimate: prop ycompnby, over(treat)

	*Substantive analysis model - risk difference
	mi estimate, mcerr vart: glm ycompnby i.treat, link(identity) family(binomial) iterate(50)
		mat def R 		= r(table)
		local b_trt 	= R[1, 2]
		local ll_trt 	= R[5, 2]
		local ul_trt 	= R[6, 2]
		local p 		= R[4, 2]

	post `topps' (1) ("MI compn by trt ice") ("Risk diff") (`b_trt') (`ll_trt') (`ul_trt') (`p')

	*Substantive analysis model - log OR
	mi estimate, mcerr vart: logit ycompnby i.treat, iterate(50)
		mat def R 		= r(table)
		local b_trt 	= R[1, 2]
		local ll_trt 	= R[5, 2]
		local ul_trt 	= R[6, 2]
		local p 		= R[4, 2]

	post `topps' (1) ("MI compn by trt ice") ("Log OR") (`b_trt') (`ll_trt') (`ul_trt') (`p')
}

else	{
    di as err "MI compn by trt ice failed with error code " as result _rc
}







/*Approach 2: 
- Each block is set to missing if bleeding status is missing for at least 3 days
- Each block takes value 1 if bleeding occurs on at least 1 day
*/

use "`outputpath'\topps_data_formi.dta", clear

local i = 1
local j = 5

forval b = 1/6	{
	tempvar block`b'_miss
	egen block`b' = rowmax(whograde`i'-whograde`j')
		
	*Check number of days with missing bleeding status within each block (overall and by treatment)
	egen `block`b'_miss' = rowmiss(whograde`i'-whograde`j')
	tab `block`b'_miss' treat, col
		
	replace block`b' = . if `block`b'_miss' > 2
		
	*Check number of participants with at least 1 bleeding episode in each block (overall and by treatment)
	tab block`b' treat, m col
		
	local i = `i' + 5
	local j = `j' + 5
}

misstable patterns block1-block6, freq

//Complete record analysis - composite y with 6 components
tempvar cra_miss
egen `cra_miss' = rowmiss(block1-block6)
gen ycra = 1 if inlist(1, block1, block2, block3, block4, block5, block6) & `cra_miss' == 0
replace ycra = 0 if missing(ycra) & `cra_miss' == 0

*Check % missing values in the composite and % of events by treatment
tempvar r_ycra
mark `r_ycra' if !missing(ycra)
tab `r_ycra' treat, col
tab ycra treat, col

*Substantive analysis model - risk difference
glm ycra i.treat, link(identity) family(binomial) iterate(50)
	mat def R 		= r(table)
	local b_trt 	= R[1, 2]
	local ll_trt 	= R[5, 2]
	local ul_trt 	= R[6, 2]
	local p 		= R[4, 2]

post `topps' (2) ("CRA") ("Risk diff") (`b_trt') (`ll_trt') (`ul_trt') (`p')

*Substantive analysis model - log OR
logit ycra i.treat, iterate(50)
	mat def R 		= r(table)
	local b_trt 	= R[1, 2]
	local ll_trt 	= R[5, 2]
	local ul_trt 	= R[6, 2]
	local p 		= R[4, 2]

post `topps' (2) ("CRA") ("Log OR") (`b_trt') (`ll_trt') (`ul_trt') (`p')

//Derived composite
gen yderiv = ycra
replace yderiv = 1 if inlist(1, block1, block2, block3, block4, block5, block6) & missing(yderiv)

*Check % missing values in the composite and % of events by treatment
tempvar r_yderiv
mark `r_yderiv' if !missing(yderiv)
tab `r_yderiv' treat, col
tab yderiv treat, col

*Substative analysis model - risk difference
glm yderiv i.treat, link(identity) family(binomial) iterate(50)
	mat def R 		= r(table)
	local b_trt 	= R[1, 2]
	local ll_trt 	= R[5, 2]
	local ul_trt 	= R[6, 2]
	local p 		= R[4, 2]

post `topps' (2) ("Deriv") ("Risk diff") (`b_trt') (`ll_trt') (`ul_trt') (`p')

*Substantive analysis model - log OR
logit yderiv i.treat, iterate(50)
	mat def R 		= r(table)
	local b_trt 	= R[1, 2]
	local ll_trt 	= R[5, 2]
	local ul_trt 	= R[6, 2]
	local p 		= R[4, 2]

post `topps' (2) ("Deriv") ("Log OR") (`b_trt') (`ll_trt') (`ul_trt') (`p')

//MI at composite level (MI-CRA)
post `postseed' (2) ("MI-CRA") (substr(c(rngstate),1,2000)) (substr(c(rngstate),2001,2000)) (substr(c(rngstate),4001,.))

local m = 50
local c = 20
compress

mi set wide
mi register imputed ycra 
mi register regular treat 

mi impute logit ycra i.treat, add(`m') augment iterate(50) dots 

*Check imputed values in the first 5 imputed datasets
forval i = 1/5	{
	tab _`i'_ycra _mi_m, col
}

*Check % of events in the composite by treatment 
mi estimate: prop ycra, over(treat)

*Substantive analysis model - risk difference
mi estimate, mcerr vart: glm ycra i.treat, link(identity) family(binomial) iterate(50)
	mat def R 		= r(table)
	local b_trt 	= R[1, 2]
	local ll_trt 	= R[5, 2]
	local ul_trt 	= R[6, 2]
	local p 		= R[4, 2]

post `topps' (2) ("MI-CRA") ("Risk diff") (`b_trt') (`ll_trt') (`ul_trt') (`p')

*Substantive analysis model - log OR
mi estimate, mcerr vart: logit ycra i.treat, iterate(50)
	mat def R 		= r(table)
	local b_trt 	= R[1, 2]
	local ll_trt 	= R[5, 2]
	local ul_trt 	= R[6, 2]
	local p 		= R[4, 2]

post `topps' (2) ("MI-CRA") ("Log OR") (`b_trt') (`ll_trt') (`ul_trt') (`p')

//MI at composite level (MI-Deriv)
post `postseed' (2) ("MI-Deriv") (substr(c(rngstate),1,2000)) (substr(c(rngstate),2001,2000)) (substr(c(rngstate),4001,.))

mi extract 0, clear
mi set wide
mi register imputed yderiv 
mi register regular treat 

mi impute logit yderiv i.treat, add(`m') augment iterate(50) dots

*Check imputed values in the first 5 imputed datasets
forval i = 1/5	{
	tab _`i'_yderiv _mi_m, col
}

*Check % of events in the composite by treatment 
mi estimate: prop yderiv, over(treat)

*Substantive analysis model - risk difference
mi estimate, mcerr vart: glm yderiv i.treat, link(identity) family(binomial) iterate(50)
	mat def R 		= r(table)
	local b_trt 	= R[1, 2]
	local ll_trt 	= R[5, 2]
	local ul_trt 	= R[6, 2]
	local p 		= R[4, 2]

post `topps' (2) ("MI-Deriv") ("Risk diff") (`b_trt') (`ll_trt') (`ul_trt') (`p')

*Substantive analysis model - log OR
mi estimate, mcerr vart: logit yderiv i.treat, iterate(50)
	mat def R 		= r(table)
	local b_trt 	= R[1, 2]
	local ll_trt 	= R[5, 2]
	local ul_trt 	= R[6, 2]
	local p 		= R[4, 2]

post `topps' (2) ("MI-Deriv") ("Log OR") (`b_trt') (`ll_trt') (`ul_trt') (`p')

//MI at component level (MI compn main)
post `postseed' (2) ("MI compn main") (substr(c(rngstate),1,2000)) (substr(c(rngstate),2001,2000)) (substr(c(rngstate),4001,.))

mi extract 0, clear
mi set wide
mi register imputed block1-block6
mi register regular treat 

mi impute chained (logit, augment iterate(50)) block1 block2 block3 block4 block5 block6 = i.treat, add(`m') burnin(`c') dots 

*Check imputed values in the e.g. 5th imputed dataset
forval i = 1/6	{
	tab _5_block`i' _mi_m, col
}

*Check % of events in each block by treatment
mi estimate: prop block1-block6, over(treat)

*Generate composite
mi passive: gen ycompnmain = inlist(1, block1, block2, block3, block4, block5, block6)

*Check passively imputed values in the first 5 imputed datasets
forval i = 1/5	{
	tab _`i'_ycompnmain _mi_m, col
}

*Check % of events in the composite by treatment 
mi estimate: prop ycompnmain, over(treat)

*Substantive analysis model - risk difference
mi estimate, mcerr vart: glm ycompnmain i.treat, link(identity) family(binomial) iterate(50)
	mat def R 		= r(table)
	local b_trt 	= R[1, 2]
	local ll_trt 	= R[5, 2]
	local ul_trt 	= R[6, 2]
	local p 		= R[4, 2]

post `topps' (2) ("MI compn main") ("Risk diff") (`b_trt') (`ll_trt') (`ul_trt') (`p')

*Substantive analysis model - log OR
mi estimate, mcerr vart: logit ycompnmain i.treat, iterate(50)
	mat def R 		= r(table)
	local b_trt 	= R[1, 2]
	local ll_trt 	= R[5, 2]
	local ul_trt 	= R[6, 2]
	local p 		= R[4, 2]

post `topps' (2) ("MI compn main") ("Log OR") (`b_trt') (`ll_trt') (`ul_trt') (`p')

cap drop ycompnmain
mi update

//MI at component level (MI compn by trt)

*MI using -mi impute chained-, conditioning on all other blocks when imputing each block (failed due to non-convergence in past runs)
post `postseed' (2) ("MI compn by trt") (substr(c(rngstate),1,2000)) (substr(c(rngstate),2001,2000)) (substr(c(rngstate),4001,.))

mi extract 0, clear
mi set wide
mi register imputed block1-block6
mi register regular treat 

cap mi impute chained (logit, augment iterate(50)) block1 block2 block3 block4 block5 block6, by(treat) add(`m') burnin(`c') 

if !_rc	{		  
	*Check imputed values in the e.g. 5th imputed dataset
	forval i = 1/6	{
		tab _5_block`i' _mi_m, col
	}

	*Check % of events in each block by treatment
	mi estimate: prop block1-block6, over(treat)

	*Generate composite
	mi passive: gen ycompnby = inlist(1, block1, block2, block3, block4, block5, block6)

	*Check passively imputed values in the first 5 imputed datasets
	forval i = 1/5	{
		tab _`i'_ycompnby _mi_m, col
	}

	*Check % of events in the composite by treatment 
	mi estimate: prop ycompnby, over(treat)

	*Substantive analysis model - risk difference
	mi estimate, mcerr vart: glm ycompnby i.treat, link(identity) family(binomial) iterate(50)
		mat def R 		= r(table)
		local b_trt 	= R[1, 2]
		local ll_trt 	= R[5, 2]
		local ul_trt 	= R[6, 2]
		local p 		= R[4, 2]

	post `topps' (2) ("MI compn by trt") ("Risk diff") (`b_trt') (`ll_trt') (`ul_trt') (`p')

	*Substantive analysis model - log OR
	mi estimate, mcerr vart: logit ycompnby i.treat, iterate(50)
		mat def R 		= r(table)
		local b_trt 	= R[1, 2]
		local ll_trt 	= R[5, 2]
		local ul_trt 	= R[6, 2]
		local p 		= R[4, 2]

	post `topps' (2) ("MI compn by trt") ("Log OR") (`b_trt') (`ll_trt') (`ul_trt') (`p')
	
	cap drop ycompnby
	mi update
}

else	{
	di as err "MI compn by trt failed with error code " as result _rc
}

*MI using -mi impute chained-, conditioning on 2 neighbouring blocks when imputing each block
post `postseed' (2) ("MI compn by trt cond2") (substr(c(rngstate),1,2000)) (substr(c(rngstate),2001,2000)) (substr(c(rngstate),4001,.))

mi extract 0, clear
mi set wide
mi register imputed block1-block6
mi register regular treat

cap mi impute chained (logit, omit(i.block4 i.block5 i.block6) augment iterate(50)) block1 ///
					  (logit, omit(i.block4 i.block5 i.block6) augment iterate(50)) block2 ///
				      (logit, omit(i.block1 i.block5 i.block6) augment iterate(50)) block3 ///
					  (logit, omit(i.block1 i.block2 i.block6) augment iterate(50)) block4 ///
				      (logit, omit(i.block1 i.block2 i.block3) augment iterate(50)) block5 ///
					  (logit, omit(i.block1 i.block2 i.block3) augment iterate(50)) block6 ///
				      , by(treat) add(`m') burnin(`c')
				  
if !_rc	{		  
	*Check imputed values in the e.g. 5th imputed dataset
	forval i = 1/6	{
		tab _5_block`i' _mi_m, col
	}

	*Check % of events in each block by treatment
	mi estimate: prop block1-block6, over(treat)

	*Generate composite
	mi passive: gen ycompnby = inlist(1, block1, block2, block3, block4, block5, block6)

	*Check passively imputed values in the first 5 imputed datasets
	forval i = 1/5	{
		tab _`i'_ycompnby _mi_m, col
	}

	*Check % of events in the composite by treatment 
	mi estimate: prop ycompnby, over(treat)

	*Substantive analysis model - risk difference
	mi estimate, mcerr vart: glm ycompnby i.treat, link(identity) family(binomial) iterate(50)
		mat def R 		= r(table)
		local b_trt 	= R[1, 2]
		local ll_trt 	= R[5, 2]
		local ul_trt 	= R[6, 2]
		local p 		= R[4, 2]

	post `topps' (2) ("MI compn by trt cond2") ("Risk diff") (`b_trt') (`ll_trt') (`ul_trt') (`p')

	*Substantive analysis model - log OR
	mi estimate, mcerr vart: logit ycompnby i.treat, iterate(50)
		mat def R 		= r(table)
		local b_trt 	= R[1, 2]
		local ll_trt 	= R[5, 2]
		local ul_trt 	= R[6, 2]
		local p 		= R[4, 2]

	post `topps' (2) ("MI compn by trt cond2") ("Log OR") (`b_trt') (`ll_trt') (`ul_trt') (`p')
	
	cap drop ycompnby
	mi update
}

else	{
	di as err "MI compn by trt cond2 failed with error code " as result _rc
}

*MI using -ice-, conditioning on all other blocks when imputing each block (uvis automatically checks for perfect prediction)
post `postseed' (2) ("MI compn by trt ice") (substr(c(rngstate),1,2000)) (substr(c(rngstate),2001,2000)) (substr(c(rngstate),4001,.))

mi extract 0, clear

cap ice block1 block2 block3 block4 block5 block6, by(treat) m(`m') cycles(`c') saving("`outputpath'\mi_compn_by_trt_ice_approach2.dta", replace) 

if !_rc	{
	use "`outputpath'\mi_compn_by_trt_ice_approach2.dta", clear

	mi import ice, imputed(block1-block6) clear
	mi convert wide, clear
													 
	*Check imputed values in the e.g. 5th imputed dataset
	forval i = 1/6	{
		tab _5_block`i' _mi_m, col
	}

	*Check % of events in each block by treatment
	mi estimate: prop block1-block6, over(treat)

	*Generate composite
	mi passive: gen ycompnby = inlist(1, block1, block2, block3, block4, block5, block6)

	*Check passively imputed values in the first 5 imputed datasets
	forval i = 1/5	{
		tab _`i'_ycompnby _mi_m, col
	}

	*Check % of events in the composite by treatment 
	mi estimate: prop ycompnby, over(treat)

	*Substantive analysis model - risk difference
	mi estimate, mcerr vart: glm ycompnby i.treat, link(identity) family(binomial) iterate(50)
		mat def R 		= r(table)
		local b_trt 	= R[1, 2]
		local ll_trt 	= R[5, 2]
		local ul_trt 	= R[6, 2]
		local p 		= R[4, 2]

	post `topps' (2) ("MI compn by trt ice") ("Risk diff") (`b_trt') (`ll_trt') (`ul_trt') (`p')

	*Substantive analysis model - log OR
	mi estimate, mcerr vart: logit ycompnby i.treat, iterate(50)
		mat def R 		= r(table)
		local b_trt 	= R[1, 2]
		local ll_trt 	= R[5, 2]
		local ul_trt 	= R[6, 2]
		local p 		= R[4, 2]

	post `topps' (2) ("MI compn by trt ice") ("Log OR") (`b_trt') (`ll_trt') (`ul_trt') (`p')
}

else	{
    di as err "MI compn by trt ice failed with error code " as result _rc
}

postclose `topps'

postclose `postseed'


log close
