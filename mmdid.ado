/* This program is to illustrate and replicate for the applications in paper
First version 22/11/2024
Please cite as: Quach and Wichitaksorn (2024) "Multilevel modeling and Difference-in-Differences"
*/

capture pro drop _all
program  mmdid, eclass byable(recall)
	syntax  varlist [if] [in] [pw fw] , /*
				*/ TVar(varlist min=1 max=1) /*	
				*/ Evtime(varlist min=1 max=1) [ /*
				*/ ATtribute(varlist numeric) /* 
				*/ TObit  /*
				*/ NOLOG  /*
				*/ LL(numlist min=1 max=1) /*
				*/ UL(numlist min=1 max=1) /*
				*/ HEADer FETable ]
	
	if "`tobit'" != "" {
		local tobit "`tobit'"
	}
	
	if "`nolog'" != "" {
		local nolog "`nolog'"
	}
	
 	if "`header'" == "" {
		local header "nohead"
	}
	else {
		local header ""
	}
	
	if "`fetable'" =="" {
		local fetable "nofet"
	}
	else {
		local fetable ""
	}
	
	qui sum `evtime'
	if r(min) >=0 {
		di as error "Please ensure that the calendar time has been converted into" //
	di as error	"event time (zero normalization), that the untreated units have the same time," //
	di as error	"and that there is a difference from the treated units."
	exit(198)
	}
	
	qui sum `tvar'
	if r(max) - r(min) > 1 {
		di as error "Please ensure that the treatment variable is defined as 0 and 1"
		exit(198)
	}
	
	
	tokenize `attribute'
	local het1  "`1'"
	local het2  "`2'"
	
	if "`attribute'" == "" {
		local count = 0
	}
	else {
		local count : word count `attribute'
	}

	qui {
		reg `varlist'
		local lhs "`e(depvar)'"
		_getrhs varlist
		local nrhs: list sizeof rhs
		local rhs "`varlist'"
		
		tempvar wt touse
		if "`weight'" != "" {
			qui gen double `wt' `exp'
			local wtopt "[`weight' = `wt']"
			local wtopt1 "`weight'(`wt')"
		} 
		else {
			qui gen byte `wt' = 1
		}
	
		mark `touse' `if' `in' `wtopt'
		markout `touse' `varlist'
		
		tempvar n nj njj njjj
		gen `n' = 1
		
		if "`attribute'"  == "" {
			egen `nj' = sum(`n') if `touse', by(`evtime')
		}
		else if `count' ==  1 {
			egen `nj' = sum(`n') if `touse' , by(`evtime')
			egen `njj' = sum(`n') if `touse' , by(`evtime' "`het1'")
		}
		else {
			egen `nj' = sum(`n') if `touse', by(`evtime')
			egen `njj' = sum(`n') if `touse', by(`evtime' "`het1'")
			egen `njjj' = sum(`n') if `touse', by(`evtime' "`het1' `het2'")
		}
		
		// Estimate
		
		if "`tobit'" == "" & "`attribute'" == "" {
		noi	xtmixed `lhs' `tvar' `rhs' `wtopt' if `touse' || `evtime':  ,	`nolog'	`wtopt1' `header' `fetable'
		}
		else if "`tobit'" != "" & "`attribute'" =="" {
		noi	metobit `lhs' `tvar' `rhs' `wtopt' if `touse' || `evtime':  ,	`nolog' `wtopt1' `header' `fetable'
		}
		else if "`tobit'" == "" & `count' ==1 {
		noi	xtmixed `lhs' `tvar' `rhs' `wtopt' if `touse' || `evtime':  || `het1': , `nolog' `wtopt1' `header' `fetable'
		}
		else if "`tobit'" != "" & `count' ==1 {
		noi	metobit `lhs' `tvar' `rhs' `wtopt' if `touse' || `evtime':  || `het1': , `nolog' `wtopt1' `header' `fetable'
		}
		else if "`tobit'" == "" & `count' ==2 {
		noi xtmixed `lhs' `tvar' `rhs' `wtopt' if `touse' || `evtime': || `het1': || `het2': , `nolog'  `wtopt1' `header' `fetable'
		}
		else {
		noi metobit `lhs' `tvar' `rhs' `wtopt' if `touse' || `evtime': || `het1': || `het2': , `nolog'  `wtopt1' `header' `fetable' ll(`ll') ul(`ul') 
		}
		
		// Predict and adjust parameters for computing ATT
		
		tempvar u1 u2 u3 ue1 ue2 ue3 u_a1 u_a2 u_a3 uh uh2
		tempname sigma_u1 sigma_u2 sigma_u3 sigma_e  pchi2
		scalar `pchi2' = e(p_c)
		if "`tobit'" == "" & "`attribute'" == "" {
			predict `u1' if `touse' , reffect
			predict `ue1' if `touse'  , reses
			gen `u_a1' = `u1' / (exp(_b[lns1_1_1:_cons])^2 	///
								/ (exp(_b[lns1_1_1:_cons])^2 + exp(_b[lnsig_e:_cons])^2 / `nj'))	if `touse' 
			
		} 
		else if "`tobit'" != "" & "`attribute'" =="" {
			predict `u1' if `touse'  , reffect reses(`ue1')
			gen `u_a1' = `u1' / (_b[/:var(_cons[`evtime'])]  ///
								/ (_b[/:var(_cons[`evtime'])] + _b[/:var(e.`lhs')] / `nj')) if `touse' 
		}
		
		else if  "`tobit'" == "" & `count' ==1 {
			predict `u1' `u2' if `touse' , reffect
			predict `ue1' `ue2' if `touse' , reses
			predict `uh' if `touse' , resid
			gen `u_a1' = `u1' / (exp(_b[lns1_1_1:_cons])^2 	///
								/ (exp(_b[lns1_1_1:_cons])^2 + exp(_b[lnsig_e:_cons])^2 / `nj'))	if `touse' 
			gen `u_a2' = `u2' / (exp(_b[lns2_1_1:_cons])^2 	///
								/ (exp(_b[lns2_1_1:_cons])^2 + exp(_b[lnsig_e:_cons])^2 / `njj'))	if `touse' 
		/*	replace `ue1' = sqrt(`ue1'^2 / (exp(_b[lns1_1_1:_cons])^2 	///
								/ (exp(_b[lns1_1_1:_cons])^2 + exp(_b[lnsig_e:_cons])^2 / `nj')))		
			replace `ue2' = sqrt(`ue2'^2 / (exp(_b[lns2_1_1:_cons])^2 	///
								/ (exp(_b[lns2_1_1:_cons])^2 + exp(_b[lnsig_e:_cons])^2 / `njj')))
			gen `uh2' = `uh'^2 */
		}
		else if "`tobit'" != "" & `count' ==1 {
			predict `u1' `u2' if `touse' , reffect reses(`ue1' `ue2')
			gen `u_a1' = `u1' / (_b[/:var(_cons[`evtime'])]  ///
								/ (_b[/:var(_cons[`evtime'])] + _b[/:var(e.`lhs')] / `nj')) if `touse' 
			gen `u_a2' = `u2' / (_b[/:var(_cons[`het1'])]  ///
								/ (_b[/:var(_cons[`het1'])] + _b[/:var(e.`lhs')] / `njj')) if `touse' 
		}
		
		else if "`tobit'" == "" & `count' ==2 {
			predict `u1' `u2' `u3' if `touse' , reffect
			predict `ue1' `ue2' `ue3' if `touse' , reses
			gen `u_a1' = `u1' / (exp(_b[lns1_1_1:_cons])^2 	///
								/ (exp(_b[lns1_1_1:_cons])^2 + exp(_b[lnsig_e:_cons])^2 / `nj')) if `touse' 
			gen `u_a2' = `u2' / (exp(_b[lns2_1_1:_cons])^2 	///
								/ (exp(_b[lns2_1_1:_cons])^2 + exp(_b[lnsig_e:_cons])^2 / `njj')) if `touse' 
			gen `u_a3' = `u3' / (exp(_b[lns3_1_1:_cons])^2 	///
								/ (exp(_b[lns3_1_1:_cons])^2 + exp(_b[lnsig_e:_cons])^2 / `njjj')) if `touse' 
		}
		else {
			predict `u1' `u2' `u3' if `touse' , reffect reses(`ue1' `ue2' `ue3')
			gen `u_a1' = `u1' /* (_b[/:var(_cons[`evtime'])]  ///
								/ (_b[/:var(_cons[`evtime'])] + _b[/:var(e.`lhs')] / `nj')) */
			gen `u_a2' = `u2' /* (_b[/:var(_cons[`evtime'>`het1'])]  ///
								/ (_b[/:var(_cons[`evtime'>`het1'])] + _b[/:var(e.`lhs')] / `njj')) */
			gen `u_a3' = `u3' /* (_b[/:var(_cons[`evtime'>`het1'>`het2'])]  ///
								/ (_b[/:var(_cons[`evtime'>`het1'>`het2'])] + _b[/:var(e.`lhs')] / `njjj')) */
		}
		
		/// Adjust standard error at level 2 for parallel trend test	
		
		if "`attribute'" =="" {
			scalar `sigma_u1' = (exp(_b[lns1_1_1:_cons]))^2	
		}
		if "`attribute'" != "" & "`tobit'" == "" {
			scalar `sigma_u1' = (exp(_b[lns1_1_1:_cons]))^2	
			scalar `sigma_u2' = (exp(_b[lns2_1_1:_cons]))^2
		}
		/*
		if "`tobit'" == "" {
			if "`attribute'"=="" &  `sigma_u1' < 1e-05 {			
				replace `ue1' = `ue1'^2 / (exp(_b[lns1_1_1:_cons])^2 	///
										/ (exp(_b[lns1_1_1:_cons])^2 + exp(_b[lnsig_e:_cons])^2 / `nj'))	
					}
				
			if `count'==1   {
					replace `ue1' = `ue1'^2 / (exp(_b[lns1_1_1:_cons])^2 	///
											/ (exp(_b[lns1_1_1:_cons])^2 + exp(_b[lnsig_e:_cons])^2 / `nj'))
					replace `ue2' = `ue2'^2 / (exp(_b[lns2_1_1:_cons])^2 	///
											/ (exp(_b[lns2_1_1:_cons])^2 + exp(_b[lnsig_e:_cons])^2 / `njj'))
			}
		}
		
		if `count'==2 {
				replace `ue1' = `ue1'^2 / (exp(_b[lns1_1_1:_cons])^2 	///
											/ (exp(_b[lns1_1_1:_cons])^2 + exp(_b[lnsig_e:_cons])^2 / `nj'))
				replace `ue2' = `ue2'^2 / (exp(_b[lns2_1_1:_cons])^2 	///
											/ (exp(_b[lns2_1_1:_cons])^2 + exp(_b[lnsig_e:_cons])^2 / `njj'))
				replace `ue3' = `ue3'^2 / (exp(_b[lns3_1_1:_cons])^2 	///
											/ (exp(_b[lns3_1_1:_cons])^2 + exp(_b[lnsig_e:_cons])^2 / `njjj'))
		}
		*/
		// Identify groups time
			
		tempvar tag1 tag2 tag3 tag22 tag33 ua ue pp ppe att atte ate atee nn 
		tempname  mhet0 mhet1 mhet2 obhet0 obhet1 obhet2  ev_min s2
	
		sum `evtime' if `touse' 
		scalar `ev_min' = r(min)
		egen `tag1' = tag(`evtime') if `touse' 
		egen `tag2' = tag(`evtime' `het1') if `touse' 
		egen `tag3' = tag(`evtime' `het1' `het2') if `touse' 
	//	sum `uh2'
		scalar `s2' =r(mean)
		if `count' >0 {
			egen `tag22' = tag(`het1' `tvar') if `touse' 
			egen `tag33' = tag(`het2' `tvar') if `touse' 
		}
		
		if `count' ==2 {
			gen `ua' = `u_a1' + `u_a2' + `u_a3' if `touse' 
			gen `ue' = `ue3'^2  + `ue2'^2 + `ue1'^2 if `touse' 
		}
		else if `count'==1 {
			gen `ua' = `u_a1' + `u_a2' if `touse' 
				if `pchi2' >0.1 {
					replace `ue1' = sqrt(`ue1'^2 / (exp(_b[lns1_1_1:_cons])^2 	///
											/ (exp(_b[lns1_1_1:_cons])^2 + exp(_b[lnsig_e:_cons])^2 / `nj'))) if `touse' 
					replace `ue2' = sqrt(`ue2'^2 / (exp(_b[lns2_1_1:_cons])^2 	///
											/ (exp(_b[lns2_1_1:_cons])^2 + exp(_b[lnsig_e:_cons])^2 / `njj'))) if `touse' 					
				}
			gen `ue' = (`ue2'^2 + `ue1'^2) if `touse' 
		}
		else {
			gen `ua' = `u_a1' if `touse' 
				if `pchi2' > 0.1 {
					replace `ue1' = sqrt( `ue1'^2 / (exp(_b[lns1_1_1:_cons])^2 	///
											/ (exp(_b[lns1_1_1:_cons])^2 + exp(_b[lnsig_e:_cons])^2 / `nj'))) if `touse' 
				}
			gen `ue' = `ue1'^2 if `touse' 
		}
		
		sort `het2' `het1' `evtime'
		
		if `count'==2 {
			mkmat `het1' if `tag22' ==1 & `tvar'==1 & `touse' 
			mat `mhet1' = `het1''
			mkmat `het2' if `tag33' ==1 & `tvar'==1 & `touse' 
			mat `mhet2' = `het2''
			sum `het1' if `tag22' == 1 & `tvar'==1 & `touse' 
			local obhet1 = r(N)
			sum `het2' if `tag33' == 1 & `tvar'==1 & `touse' 
			local obhet2 = r(N)
		}
		else if `count' ==1 {
			mkmat `het1' if `tag22' ==1 & `tvar'==1 & `touse' 
			mat `mhet1' = `het1''
			sum `het1' if `tag22' == 1 & `tvar'==1 & `touse' 
			local obhet1 = r(N)
		}
		else {
			mkmat `evtime' if `tag1' ==1 & `tvar'==1 & `touse' 
			mat `mhet0' = `evtime''
		}
	
		
		
		local numlist1 ""
		local numlist2 ""
		
		if `count' ==2 {
			forvalues i = 1/`obhet1' {
				
				local numlist1 `numlist1' `=`mhet1'[1,`i']'
				}
			forvalues i = 1/`obhet2' {
			
				local numlist2 `numlist2' `=`mhet2'[1,`i']'
				}
		}
		else if `count'==1 {
			forvalues i = 1/`obhet1' {
			//	local numlist1 "`numlist1' `" "'"
				local numlist1 `numlist1' `=`mhet1'[1,`i']'
			}
		}
		
		
		/// Calculate counterfactual
		gen `pp' =.
		gen `ppe' =.
	
		
		if `count' ==2 {
			foreach i in `numlist1' {
				foreach j in `numlist2' {
					sum `ua' if `evtime' <0 & `tvar'==1 & `het1' ==`i' & `het2' ==`j' & `touse' 
					replace `pp' = r(mean) if `het1'==`i' & `het2'==`j' & `touse' 
					sum `ue' if `evtime' <0 & `tvar'==1 & `het1' ==`i' & `het2' ==`j' & `touse' 
					replace `ppe' =r(mean) / r(N) if `het1'==`i' & `het2'==`j' & `touse' 
				} 
			}
		}
		else if `count' ==1 {
			foreach i of local numlist1 {
				sum `ua' if `evtime' <0 & `tvar'==1 & `het1' ==`i' & `touse' 
				replace `pp' = r(mean) if `het1'==`i' & `touse' 
				sum `ue' if `evtime' <0 & `tvar'==1 & `het1' ==`i' & `touse'  
				replace `ppe' = r(mean) / r(N) if  `het1'==`i' & `touse' 
			}
		}
		else {
			sum `ua' if `evtime' <0 & `tvar'==1 & `touse' 
			replace `pp' = r(mean) if `touse' 
			sum `ue' if  `evtime' <0 & `tvar'==1 & `touse' 
			replace `ppe' = r(mean) / r(N) if `touse' 
		}
		
		
		/// Get name of column for displaying
		
	
		local colname 
		
		if `count' ==1 {
		
		sort `het1' `evtime' 
		foreach i of local numlist1 {
			sum `evtime' if `het1'==`i' & `tvar'==1 & `touse'
				tempname rmin rmax
				local rmin =r(min)
				local rmax =r(max)
				forv j = `rmin'(1)`rmax' {
					local colname `colname' `i':`j'
				}
			}
		}
		else if `count' ==2 {
			foreach c of local numlist2 {
				sort  `het1' `evtime'
				foreach i of local numlist1 {
					sum `evtime' if `het1'==`i' & `tvar'==1 & `touse'
					tempname rmin rmax
					local rmin =r(min)
					local rmax =r(max)
						forv j = `rmin'(1)`rmax' {
							local colname `colname' `c':`i':`j'
					}
				}
			}
		}
		else {
			sum `evtime' if `tvar'==1 & `touse'
			tempname rmin rmax
			local rmin = r(min)
			local rmax = r(max)
				forv i = `rmin'(1)`rmax' {
					local colname `colname' Event_time:`i'
				}
		}
		
		
		// Calculate attribute-cluster ATT
		tempvar atte1
		gen `att' = `ua' - `pp' if `touse' 
		label var `att' "ATT"
		gen `atte' = `ue' + `ppe' if `touse' 
		gen `atte1' = sqrt(`atte') if `touse' 
	//	tw scatter `att' `evtime' if `tag3'==1, by(`het1' )
		
		tempname result stdre ei attc Ei Attc Result Std Attc2 tag111
		
		sort `het2' `het1' `evtime'
	
		if `count' > 0 {
			mkmat `evtime' if `tag3' ==1 & `tvar'==1 & `touse' , mat(`ei')
			mkmat `het1' if `tag3'==1 &  `tvar'==1 & `touse', mat(`attc')
			mkmat `att' if `tag3'==1 &  `tvar'==1 & `touse', mat(`result')
			mkmat `atte' if `tag3'==1 &  `tvar'==1 & `touse', mat(`stdre')
			
			mkmat `evtime' if `tag3' ==1 & `tvar'==1 & `touse' , mat(`Ei')
			mkmat `het1' if `tag3'==1 &  `tvar'==1 & `touse', mat(`Attc')
			mkmat `att' if `tag3'==1 &  `tvar'==1 & `touse', mat(`Result')
			mkmat `atte' if `tag3'==1 &  `tvar'==1 & `touse', mat(`Std')
			if `count' == 2 {
				mkmat `het2' if `tag3'==1 &  `tvar'==1 & `touse', mat(`Attc2')
			}
		}
		else {
			mkmat `evtime' if `tag3' ==1 & `tvar'==1 & `touse' , mat(`ei')
			mkmat `att' if `tag3'==1 & `tvar'==1 & `touse', mat(`result')
			mkmat `atte' if `tag3'==1 & `tvar'==1 & `touse' , mat(`stdre')
			
			mkmat `evtime' if `tag3' ==1 & `tvar'==1 & `touse' , mat(`Ei')
			mkmat `att' if `tag3'==1 & `tvar'==1 & `touse', mat(`Result')
			mkmat `atte' if `tag3'==1 & `tvar'==1 & `touse' , mat(`Std')
			
	
			
		}
	
		// Calculate final ATT 
	/*	
		local k = 1
		local uq
		 {
		foreach i of local numlist1 {
		
			sum `evtime' if `het1'==`i' & `tvar'==1 & `touse'
					tempname rmin rmax
					local rmin =r(min)
					local rmax =r(max)
					forv j = `rmin'(1)`rmax' {
						tempvar u_`i'_`k'
						gen `u_`i'_`k'' = cond(`het1'==`i' & `evtime'==`j' & `tvar'==1, 1,0)
						local uq `uq' `u_`i'_`k''
						local ++k
					}
		}
		tempname M 
		mat accum `M' =  `uq' [iweight=`uh2'] if `touse', nocons
		}
	
		*/
		
		mat `result' = `result''

	
		mat `stdre' = diag(`stdre') 
	
	//	mat `stdre' = (1/`s2')*`stdre'
	//	mat `stdre' = `stdre'*`M'*`stdre' 
		
		matrix colname `result' = `colname'
		matrix colname `stdre' = `colname'
		matrix rownames `stdre' = `colname'
	/*	eret post `result' `stdre'
		di ""
		if `count'==2 {
			di "Attribute-cluster ATT by `het1' and  `het2' "
		}
		else if `count'==1 {
			di "Attribute-cluster ATT by  `het1' "
		}
		eret display */
		
		tempvar nn3 ateh hnn hn atteh attehh tag32
		tempname result_h stdre_h eih ei_hh result_hh stdre_hh hh
		
		if `count' > 1 {
		
		egen `ateh' = mean(`att') if `tvar'==1 &`touse' , by(`evtime' `het2')
		egen `hn' = sum(`n') if `tvar'==1 & `touse', by(`evtime' `het2')
		egen `hnn' = sum(`n') if `tvar'==1 & `touse', by(`evtime' `het1')
		gen `atteh' = `atte' * `hnn' / `hn' if `touse'
		egen `attehh' = mean(`atteh') if `tvar'==1 & `touse', by(`evtime' `het2')
	
		sort `het2'  `evtime'
		egen `tag32' = tag(`evtime' `het2') if `touse'
		mkmat `evtime' if `tag32'==1 & `tvar'==1 & `touse' , mat(`eih')
		mkmat `ateh' if `tag32'==1 & `tvar'==1 & `touse', mat(`result_h')
		mkmat `attehh' if `tag32'==1 & `tvar'==1 & `touse', mat(`stdre_h')
		
		mkmat `evtime' if `tag32'==1 & `tvar'==1 & `touse' , mat(`ei_hh')
		mkmat `ateh' if `tag32'==1 & `tvar'==1 & `touse', mat(`result_hh')
		mkmat `attehh' if `tag32'==1 & `tvar'==1 & `touse', mat(`stdre_hh')
		mkmat `het2' if `tag32'==1 & `tvar'==1 & `touse', mat(`hh')
		
		mat `result_h' = `result_h''
		mat `stdre_h' = diag(`stdre_h')
		local colnameh
		foreach c of local numlist2 {
			sort `het1' `evtime'
			sum `evtime' if `het2'==`c' & `tvar'==1 & `touse'
			tempname rmin rmax
					local rmin =r(min)
					local rmax =r(max)
						forv i = `rmin'(1)`rmax' {
							local colnameh `colnameh' `c':`i'
						}
		}
	
		matrix colname `result_h' = `colnameh'
		matrix colname `stdre_h' = `colnameh'
		matrix rownames `stdre_h' = `colnameh'
		}
		
		egen `ate' = mean(`att') if  `tvar'==1 & `touse' , by(`evtime')
		egen `nn' = sum(`n') if  `tvar'==1 & `touse' , by(`evtime')
		egen `nn3' = sum(`n') if `tvar'==1 & `touse' , by(`evtime' `het1' `het2')
		replace `atte' = `atte' * `nn3' / `nn' if `touse'
		egen `atee' = mean(`atte') if  `tvar'==1 & `touse', by(`evtime')
	//	replace `atee' = `atee' / `nn'
		
		tempname result1 stdre1 ei1 Ei_f Result_f Std_f
	
		if `count' > 0 {
				sort `evtime' 
				mkmat `evtime' if `tag1' == 1 & `tvar'==1 & `touse' , mat(`ei1')
				mkmat `ate' if `tag1'==1 & `tvar'==1 & `touse' , mat(`result1')
				mkmat `atee' if `tag1'==1 & `tvar'==1 & `touse' , mat(`stdre1')
				
				mkmat `evtime' if `tag1' == 1 & `tvar'==1 & `touse' , mat(`Ei_f')
				mkmat `ate' if `tag1'==1 & `tvar'==1 & `touse' , mat(`Result_f')
				mkmat `atee' if `tag1'==1 & `tvar'==1 & `touse' , mat(`Std_f')
				
				

				
		/*		local k = 1
				local uq
			sum `evtime' if  `tvar'==1 & `touse'
					tempname rmin rmax
					local rmin =r(min)
					local rmax =r(max)
					forv j = `rmin'(1)`rmax' {
						tempvar u_`i'_`k'
						gen `u_`i'_`k'' = cond( `evtime'==`j' & `tvar'==1, 1,0)
						local uq `uq' `u_`i'_`k''
						local ++k
					}
			tempname M1
			mat accum `M1' =  `uq' [iweight=`uh2'] if `touse', nocons		
		
				
			*/	
				mat `result1' = `result1'' 
				mat `stdre1' = diag(`stdre1') 
			//	mat `stdre1' = `stdre1'*`M1'*`stdre1'
				sum `evtime' if `tvar'==1 & `touse'
				tempname rmin rmax
				local rmin = r(min)
				local rmax = r(max)
				local colname1
				forv i = `rmin'(1)`rmax' {
					local colname1 `colname1' Event_time:`i'
				}
				
				matrix colname `result1' = `colname1'
				matrix colname `stdre1' = `colname1'
				matrix rownames `stdre1' = `colname1'
			//	eret post `result1' `stdre1'
			//	eret display 
		}
		
	}
		
		eret post `result' `stdre'  
		di ""
		if `count'==2 {
			di "Attribute-cluster ATT by `het1' and  `het2' "
		}
		else if `count'==1 {
			di "Attribute-cluster ATT by  `het1' "
		}
		else if `count'==0 {
			di "ATT by event time"
		}
		eret display
		
		if `count' > 0 {
			if `count'==2 {
			di "Final aggregate ATT by `het2'"
			eret post `result_h'  `stdre_h'  
			eret dis  
			}
			di "Final aggregate ATT by event time"
			eret post `result1' `stdre1' 
			eret display  
		}
		
		
//	table `het2' `evtime' `het1' , stat(mean `ua' `ue')
	 
//	 gen ua = `ua'
//	 gen ue = `ue'
		if `count' == 0 {
			eret mat ei = `Ei'	
			eret mat result = `Result'
			eret mat stdre = `Std'
		}
		if `count' > 0 {
			eret mat ei = `Ei'
			eret mat result = `Result'
			eret mat stdre = `Std'
			eret mat attc = `Attc'
			if `count' ==2 {
				eret mat attc2 = `Attc2'
				eret mat result_h = `result_hh'
				eret mat stdre_h = `stdre_hh'
				eret mat ei_h = `ei_hh'
				eret mat hh=`hh'
				local het2 "`het2'"
				eret local att_name "`het2'"
				}
			eret mat ei_f = `Ei_f'
			eret mat result_f = `Result_f'
			eret mat stdre_f = `Std_f'
			eret local treatvar "`tvar'"
			
		}
		
end


pro mmdid_plot,
syntax [, SEParate HETerogeneity ]
 qui {
	tempname evtime attribute att att_f att_se  result_f stdre_f evtime_f attribute2 result_h stdre_h ei_h hh
	tempvar  att att_se  ei att_cl ATT SE ATT_f SE_f ei_f att_cl2 ATT_h SE_h ei_h hh
	
	local att_name `e(att_name)'
	local treatvar e(treatvar)
	mat `evtime' = e(ei)
	mat `attribute' = e(attc)
	mat `attribute2' = e(attc2)
	mat `att' = e(result)
	mat `att_se' =e(stdre)
	mat `evtime_f' = e(ei_f)
	mat `result_f' = e(result_f)
	mat `stdre_f' = e(stdre_f)
	mat `result_h' = e(result_h)
	mat `stdre_h' = e(stdre_h)
	mat `ei_h' = e(ei_h)
	mat `hh' = e(hh)
	mat list `ei_h'
	svmat `evtime', name(`ei')
	svmat `attribute', name(`att_cl')
	svmat `attribute2' , name(`att_cl2')
	svmat `att' , name(`ATT')
	svmat `att_se' , name(`SE')
	svmat `evtime_f', name(`ei_f')
	svmat `result_f', name(`ATT_f')
	svmat `stdre_f' , name(`SE_f')
	svmat `result_h', name(`ATT_h')
	svmat `stdre_h', name(`SE_h')
	svmat `ei_h', name(`ei_h')
	svmat `hh', name(`hh')
	
	tempvar att_ul att_ll ecol att_ul_f att_ll_f att_ul_h att_ll_h
	
	label var `att_cl' "Attribute clusters"
	
	replace `SE' = sqrt(`SE')
	replace `SE_f' = sqrt(`SE_f')
	replace `SE_h' = sqrt(`SE_h')
	
	gen `att_ul' = `ATT' + 1.96*`SE'
	gen `att_ll' = `ATT' - 1.96*`SE'
	
	gen `att_ul_f' = `ATT_f' + 1.96*`SE_f'
	gen `att_ll_f' = `ATT_f' - 1.96*`SE_f'
	
	gen `att_ll_h' = `ATT_h' - 1.96*`SE_h'
	gen `att_ul_h' = `ATT_h' + 1.96*`SE_h'
	
	if `att_cl2' != . {
		sum `att_cl2'
		tempname r1 r2
		scalar `r1' = r(min)
		scalar `r2' = r(max)
		local r1 = r(min)
		local r2 = r(max)
		replace `ei' = `ei' - 0.15 if `att_cl2' == `r1'
		replace `ei' = `ei' + 0.15 if `att_cl2' == `r2'
	}

	if "`heterogeneity'" != "" {
		replace `ei_h' = `ei_h' - 0.15 if `hh' == `r1'
		replace `ei_h' = `ei_h' + 0.15 if `hh' == `r2'
	}
}
	

	if "`separate'" == "" & "`heterogeneity'"=="" & `att_cl' == . {
		
		tw rbar `att_ul' `att_ll' `ei' , barwidth(0.2) lwidth(none)  || ///
			scatter `ATT' `ei' , yline(0, lpattern(dash)) xtitle("Event time") ///
								ytitle("ATT") legend(label(1 "95% CI") label(2 "ATT")) ///
								
	}
	else if "`separate'" == "" & "`heterogeneity'"=="" & `att_cl' != . {
		
		tw rbar `att_ul_f' `att_ll_f' `ei_f' , barwidth(0.2) lwidth(none)  || ///
			scatter `ATT_f' `ei_f' ,  yline(0, lpattern(dash)) xtitle("Event time") ///
								ytitle("ATT") legend(label(1 "95% CI") label(2 "ATT")) ///
								
	}
	 else if "`separate'" != "" &  `att_cl2' != . {
		
		tw rbar `att_ul' `att_ll' `ei' if `att_cl2' == `r1' , barwidth(0.2) lwidth(none)  || ///
			rbar `att_ul' `att_ll' `ei' if `att_cl2' == `r2' , barwidth(0.2) lwidth(none)  || ///
			scatter `ATT' `ei' if `att_cl2' ==`r1' , by(`att_cl') connect(L) mcolor(stc5) lcolor(dknavy)  ||  ///
			scatter `ATT' `ei' if `att_cl2' ==`r2' , by(`att_cl') yline(0,  lpattern(dash)) xtitle("Event time") ///
								ytitle("ATT") legend(label(1 "95% CI") label(2 "95% CI") label(3 "`att_name' = `r1'") ///
								label(4 "`att_name' = `r2'"))  ///
								mcolor(dkgreen) connect(L) lcolor(dknavy) msymbol(triangle)
	}
	
	 else if "`heterogeneity'" !="" & `att_cl2' !=.  {
		
		tw rbar `att_ul_h' `att_ll_h' `ei_h' if `hh' == `r1' , barwidth(0.2) lwidth(none)  || ///
			rbar `att_ul_h' `att_ll_h' `ei_h' if `hh' == `r2' , barwidth(0.2) lwidth(none)  || ///
			scatter `ATT_h' `ei_h' if `hh' ==`r1' ,  connect(L) mcolor(stc5) lcolor(dknavy)  ||  ///
			scatter `ATT_h' `ei_h' if `hh' ==`r2' ,  yline(0,  lpattern(dash)) xtitle("Event time") ///
								ytitle("ATT") legend(label(1 "95% CI") label(2 "95% CI") label(3 "`att_name' = `r1'") ///
								label(4 "`att_name' = `r2'"))  ///
								mcolor(dkgreen) connect(L) lcolor(dknavy) msymbol(triangle)
								
	}
	else {
		
		tw rbar `att_ul' `att_ll' `ei' , barwidth(0.2) lwidth(none)  || ///
			scatter `ATT' `ei' , by(`att_cl')  yline(0, lpattern(dash)) xtitle("Event time") ///
								ytitle("ATT") legend(label(1 "95% CI") label(2 "ATT")) connect(L) ///
								lcolor(dknavy)
	
		
	}
end














