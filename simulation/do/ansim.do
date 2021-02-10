//Outcome types: 1-simple; 2-complex
//Cases: 1-case I (noint); 2-case II (intctr); 3-case III (inttrtctr)
//Missingness mechanisms: 1-MCAR; 2-MAR1; 3-MAR2

*** ANALYSIS OF SIMULATION RESULTS ***

* Installation of -simsum- is required.
********************************************************************************

cd "\\ad.ucl.ac.uk\home9\zcakf79\DesktopSettings\Desktop\20210204\stata_output"

//Analyse datasets containing simulation results
foreach type in simple complex {
	foreach case in noint intctr inttrtctr	{

		tempfile ansim0
		tempfile ansimx
		
		//True parameter values
		local b0 = 0.3
		local bx = 1.35

		local mismech = 1
		foreach miss in mcar mar1 mar2	{
			use `type'_`case'_`miss'_2, clear
			replace rep = rep + 1000
			append using `type'_`case'_`miss'_1

			encode method, gen(meth)
			recode meth (3=1) (1=2) (2=3) (4=4) (5=5) (6=6) (7=7) (8=8)

			lab drop meth
			lab def meth 1 "Full data" 2 "CRA" 3 "Deriv" 4 "MI–CRA" 5 "MI–Deriv" ///
						 6 "MIC–main" 7 "MIC–x" 8 "MIC–x–z{subscript:1}" 
			lab val meth meth

			drop method
			compress			 
					
			foreach param in 0 x 	{   
			simsum b`param', true(`b`param'') methodvar(meth) ref("Full data") 	///
					   id(rep) se(se`param') df(df`param') modelsemethod(mean) mcse 	///
					   saving(`type'_`case'_`miss'_b`param', replace)
			}
			
			foreach param in 0 x	{
			use `type'_`case'_`miss'_b`param', clear
			
			drop perfmeascode
			forval j = 1/8	{
				rename b`param'`j' b`j'
				rename b`param'`j'_mcse mcse`j'
			}

			reshape long b mcse, i(perfmeasnum) j(method)
			keep if inlist(perfmeasnum, 3, 4, 7, 9)

			lab def method 1 "Full data" 2 "CRA" 3 "Deriv" 4 "MI–CRA" 5 "MI–Deriv" ///
						   6 "MIC–main" 7 "MIC–x" 8 "MIC–x–z{subscript:1}" 
			lab val method method

			gen lb = b + (invnormal(0.025)*mcse)
			gen ub = b + (invnormal(0.975)*mcse)
			drop mcse	
					 
			gen mismech = `mismech'

			cap append using `ansim`param''
			save `ansim`param'', replace
			}
			local ++mismech
		}

		*Save results into separate datasets for graphs
		use `ansimx', clear
		gen param = 1
		append using `ansim0'
		replace param = 0 if missing(param)

		save `type'_`case'_ansim, replace
		
		//Graph the results in terms of bias, SEs and coverage	
		foreach param in 0 x	{
			if "`param'" == "0" local p = 0
			else local p = 1
			
			* BIAS
			use `type'_`case'_ansim, clear
			keep if perfmeasnum == 3
			cap lab drop method

			drop if inlist(method, 4, 5)
			replace method = method - 2 if method > 4

			lab def method 1 "Full data" 2 "CRA" 3 "Deriv" 4 "MI–CRA" 5 "MI–Deriv" ///
						   6 "MIC–main" 7 "MIC–x" 8 "MIC–x–z{subscript:1}" 
			lab val method method

			lab def param 0 "{&beta}{sub:0}" 1 "{&beta}{sub:x}"
			lab val param param 

			lab def mismech 1 "MCAR" 2 "MAR1" 3 "MAR2" 
			lab val mismech mismech

			sort mismech method param
				
				twoway (scatter method b if method == 1 & param == `p', msym(o) msiz(large) mcol(mrcgrey)) 					///
					   (rspike ub lb method if method == 1 & param == `p', hor lcol(mrcgrey)) 								///
					   (scatter method b if method == 2 & param == `p', msym(o) msiz(large) mcol(mrcblue)) 					///
					   (rspike ub lb method if method == 2 & param == `p', hor lcol(mrcblue)) 								///
					   (scatter method b if method == 3 & param == `p', msym(o) msiz(large) mcol(mrcmagenta)) 				///
					   (rspike ub lb method if method == 3 & param == `p', hor lcol(mrcmagenta)) 							///
					   (scatter method b if method == 4 & param == `p', msym(o) msiz(large) mcol(mrcorange)) 				///
					   (rspike ub lb method if method == 4 & param == `p', hor lcol(mrcorange)) 							///
					   (scatter method b if method == 5 & param == `p', msym(o) msiz(large) mcol(mrcgreen)) 				///
					   (rspike ub lb method if method == 5 & param == `p', hor lcol(mrcgreen)) 								///
					   (scatter method b if method == 6 & param == `p', msym(o) msiz(large) mcol(mrcyellow))				///
					   (rspike ub lb method if method == 6 & param == `p', hor lcol(mrcyellow)) 							///
						, by(mismech, col(1) noiylabel note("") legend(off) tit("Bias in point estimate", size(medsmall)) 	///
						graphregion(margin(marginstyle))) ylab(0.5(1)6.5, tstyle(none)) ymtick(1(1)6, tstyle(none)) 		///
						yscale(reverse lstyle(none)) ytit("") xlab(, gmin gmax) xline(0, lcol(gs12) lw(vthin) noextend) 	///
						xtit("") subtitle(, size(medsmall)) 																///
						legend(order(1 "Full data" 3 "CRA" 5 "Deriv" 7 "MIC–main" 9 "MIC–x" 11 "MIC–x–z{subscript:1}") ///
						col(1) size(vsmall))

				graph save `type'_`case'_bias_b`param'.gph, replace

			*SE
			use `type'_`case'_ansim, clear
			keep if inlist(perfmeasnum, 4, 7)
			cap lab drop method

			drop if inlist(method, 4, 5)
			replace method = method - 2 if method > 4

			lab def method 1 "Full data" 2 "CRA" 3 "Deriv" 4 "MI–CRA" 5 "MI–Deriv" ///
						   6 "MIC–main" 7 "MIC–x" 8 "MIC–x–z{subscript:1}"  
			lab val method method

			lab def param 0 "{&beta}{sub:0}" 1 "{&beta}{sub:x}"
			lab val param param 

			lab def mismech 1 "MCAR" 2 "MAR1" 3 "MAR2" 
			lab val mismech mismech


			sort mismech method param

			li b mismech param if method == 1 & per == 4

			gen method1 = method
			replace method = method - 0.3
			replace method = method + 0.4 if perfmeasnum==7

			twoway (scatter method b if perfmeasnum == 4 & method1 == 1 & param == `p', msym(o) msiz(large) mcol(mrcgrey)) 			///
				   (rspike ub lb method if perfmeasnum == 4 & method1 == 1 & param == `p', hor lcol(mrcgrey)) 						///
				   (scatter method b if perfmeasnum == 7 & method1 == 1 & param == `p', msym(oh) msiz(large) mcol(mrcgrey)) 		///
				   (rspike ub lb method if perfmeasnum == 7 & method1 == 1 & param == `p', hor lcol(mrcgrey)) 						///
				   (scatter method b if perfmeasnum == 4 & method1 == 2 & param == `p', msym(o) msiz(large) mcol(mrcblue)) 			///
				   (rspike ub lb method if perfmeasnum == 4 & method1 == 2 & param == `p', hor lcol(mrcblue)) 						///
				   (scatter method b if perfmeasnum == 7 & method1 == 2 & param == `p', msym(oh) msiz(large) mcol(mrcblue)) 		///
				   (rspike ub lb method if perfmeasnum == 7 & method1 == 2 & param == `p', hor lcol(mrcblue)) 						///
				   (scatter method b if perfmeasnum == 4 & method1 == 3 & param == `p', msym(o) msiz(large) mcol(mrcmagenta)) 		///
				   (rspike ub lb method if perfmeasnum == 4 & method1 == 3 & param == `p', hor lcol(mrcmagenta)) 					///
				   (scatter method b if perfmeasnum == 7 & method1 == 3 & param == `p', msym(oh) msiz(large) mcol(mrcmagenta)) 		///
				   (rspike ub lb method if perfmeasnum == 7 & method1 == 3 & param == `p', hor lcol(mrcmagenta)) 					///
				   (scatter method b if perfmeasnum == 4 & method1 == 4 & param == `p', msym(o) msiz(large) mcol(mrcorange)) 		///
				   (rspike ub lb method if perfmeasnum == 4 & method1 == 4 & param == `p', hor lcol(mrcorange)) 					///
				   (scatter method b if perfmeasnum == 7 & method1 == 4 & param == `p', msym(oh) msiz(large) mcol(mrcorange)) 		///
				   (rspike ub lb method if perfmeasnum == 7 & method1 == 4 & param == `p', hor lcol(mrcorange)) 					///
				   (scatter method b if perfmeasnum == 4 & method1 == 5 & param == `p', msym(o) msiz(large) mcol(mrcgreen)) 		///
				   (rspike ub lb method if perfmeasnum == 4 & method1 == 5 & param == `p', hor lcol(mrcgreen)) 						///
				   (scatter method b if perfmeasnum == 7 & method1 == 5 & param == `p', msym(oh) msiz(large) mcol(mrcgreen)) 		///
				   (rspike ub lb method if perfmeasnum == 7 & method1 == 5 & param == `p', hor lcol(mrcgreen)) 						/// 
				   (scatter method b if perfmeasnum == 4 & method1 == 6 & param == `p', msym(o) msiz(large) mcol(mrcyellow)) 		///
				   (rspike ub lb method if perfmeasnum == 4 & method1 == 6 & param == `p', hor lcol(mrcyellow)) 					///
				   (scatter method b if perfmeasnum == 7 & method1 == 6 & param == `p', msym(oh) msiz(large) mcol(mrcyellow)) 		///
				   (rspike ub lb method if perfmeasnum == 7 & method1 == 6 & param == `p', hor lcol(mrcyellow)) 					/// 	
				   , by(mismech, col(1) noiylabel note("") legend(off) tit("Empirical and model SEs", size(medsmall)) 				///
				   graphregion(margin(marginstyle)) ) ylab(0.5(1)6.5, tstyle(none)) ymtick(1(1)6, tstyle(none)) 					///
				   yscale(reverse lstyle(none)) ytit("") xtit("") subtitle(, size(medsmall))										///
				   legend(order(1 "Full data" 3 "CRA" 5 "Deriv" 7 "MIC–main" 9 "MIC–x" 11 "MIC–x–z{subscript:1}") ///
				   col(1) size(vsmall))

			graph save `type'_`case'_se_b`param'.gph, replace
			
			
			*CV	
			use `type'_`case'_ansim, clear
			keep if perfmeasnum == 9
			cap lab drop method

			drop if inlist(method, 4, 5)
			replace method = method - 2 if method > 4

			lab def meth 1 "Full data" 2 "CRA" 3 "Deriv" 4 "MI–CRA" 5 "MI–Deriv" ///
						 6 "MIC–main" 7 "MIC–x" 8 "MIC–x–z{subscript:1}" 
			lab val method method

			lab def param 0 "{&beta}{sub:0}" 1 "{&beta}{sub:x}"
			lab val param param 

			lab def mismech 1 "MCAR" 2 "MAR1" 3 "MAR2" 
			lab val mismech mismech

			sort mismech method param

			sort mismech method
	
			twoway (scatter method b if method == 1 & param == `p', msym(o) msiz(large) mcol(mrcgrey)) 					///
				   (rspike ub lb method if method == 1 & param == `p', hor lcol(mrcgrey)) 								///
				   (scatter method b if method == 2 & param == `p', msym(o) msiz(large) mcol(mrcblue)) 					///
				   (rspike ub lb method if method == 2 & param == `p', hor lcol(mrcblue)) 								///
				   (scatter method b if method == 3 & param == `p', msym(o) msiz(large) mcol(mrcmagenta)) 				///
				   (rspike ub lb method if method == 3 & param == `p', hor lcol(mrcmagenta)) 							///
				   (scatter method b if method == 4 & param == `p', msym(o) msiz(large) mcol(mrcorange)) 				///
				   (rspike ub lb method if method == 5 & param == `p', hor lcol(mrcorange)) 							///
				   (scatter method b if method == 5 & param == `p', msym(o) msiz(large) mcol(mrcgreen)) 				///
				   (rspike ub lb method if method == 5 & param == `p', hor lcol(mrcgreen)) 								///
				   (scatter method b if method == 6 & param == `p', msym(o) msiz(large) mcol(mrcyellow)) 				///
				   (rspike ub lb method if method == 6 & param == `p', hor lcol(mrcyellow)) 							///
				   , by(mismech, col(1) noiylabel note("") legend(off) tit("Coverage of 95% CIs", size(medsmall)) 		///
				   graphregion(margin(marginstyle))) ylab(0.5(1)6.5, tstyle(none)) ymtick(1(1)6, tstyle(none)) 			///
				   yscale(reverse lstyle(none)) ytit("") xline(95, lcol(gs12) lw(vthin) noextend) xlab(0(20)100) 		///
				   xtit("") subtitle(, size(medsmall)) 																	///
				   legend(order(1 "Full data" 3 "CRA" 5 "Deriv" 7 "MIC–main" 9 "MIC–x" 11 "MIC–x–z{subscript:1}") ///
				   col(1) size(vsmall))
												
			graph save `type'_`case'_cv_b`param'.gph, replace
			
			*Create a blank graph to make space for legend on RHS of graph
			twoway scatteri 1 1,               ///
				   msymbol(i)                  ///
				   ylab("") xlab("")           ///
				   ytitle("") xtitle("")       ///
				   yscale(off) xscale(off)     ///
				   plotregion(lpattern(blank)) ///
				   name(blank, replace)	
			
			grc1leg2 `type'_`case'_bias_b`param'.gph `type'_`case'_se_b`param'.gph `type'_`case'_cv_b`param'.gph blank ///
				, col(4) legendfrom(`type'_`case'_bias_b`param'.gph) ring(0) pos(3) 
			graph save `type'_`case'_ansim_b`param'.gph, replace
			*graph export `type'_`case'_ansim_b`param'.pdf, replace
		}
	}	
}

/*
//Check the number of simulation repetitions with perfect prediction during MI
foreach case in noint intctr inttrtctr	{
    foreach miss in mcar mar1 mar2	{
	    noi di as err "`case' ; " "`miss'"
	    
		use complex_`case'_`miss'_2, clear
		replace rep = rep+2000
		append using complex_`case'_`miss'_2
		
		tab method errcode
	}
}
*/
