program define kstest2w, rclass sortpreserve
// This is a modified version of ksmirnov.ado (version 3.5.2) that allows for weighted observations
// Code different from ksmirnov is marked with "NEW". 
// Syntax: kstest2w y, by(d) [weights(w)] - y = outcome, d = treatment indicator, w = weight

// By Jianfei Cao & Tetsuya Kaji
// Last updated: Oct 29, 2023

	version 6, missing
	syntax varname [=/exp] [if] [in] [, /*
		*/ by(varname) weights(varname) Onesamp Exact ]
	marksample touse
	markout `touse' `by', strok
	
	// varname string not allowed
	confirm numeric variable `varlist'
	
	///// NEW - default weights are set to 1 /////
	if "`weights'" == "" {
		tempvar weights
		qui gen double `weights' = 1
	}

	// count returns result, should be put after -gen double `eval' = `exp'
	// so that `exp' could include r() from users' command
	qui count if `touse'
	if !r(N) {
		error 2000
	}
	
	if "`by'"!="" {
		if `"`exp'"'!="" { error 198 }
		qui tab `by' if `touse'
		if r(r) != 2 {
			di in red "`by' takes on " r(r) " values, not 2"
			exit 450
		}
		tempvar M V1 V2 W1 W2 // NEW - specify tempvar including weights
		sort `touse' `by'
		qui by `touse' `by': gen byte `M' = cond(_n==1,1,.) if `touse'
		sort `M' `by'
		qui gen double `V1' = `varlist' if `by'==`by'[1] & `touse'
		qui gen double `V2' = `varlist' if `by'==`by'[2] & `touse'
		local vn `V1'
		local exp `V2'
		
		///// NEW - Define weights /////
		qui gen double `W1' = `weights' if `by'==`by'[1] & `touse'
		qui gen double `W2' = `weights' if `by'==`by'[2] & `touse'
		local w1 `W1'
		local w2 `W2'
		
		///// NEW - Compute the effective sample sizes /////
		egen sum1 = sum(`w1')
		egen sumOfSquares1 = sum(`w1' ^ 2)
		local n1 = (sum1)^2 / sumOfSquares1
		egen sum2 = sum(`w2')
		egen sumOfSquares2 = sum(`w2' ^ 2)
		local n2 = (sum2)^2 / sumOfSquares2
		drop sum1 sum2 sumOfSquares1 sumOfSquares2
		
		local typ : type `by'
		if bsubstr("`typ'",1,3)=="str" { 
			local name1 = udsubstr(`by'[1],1,16)
			local name2 = udsubstr(`by'[2],1,16)
		}
		else {
			local junk = `by'[1]
			local jname1 : label (`by') `junk' 16
			local junk = `by'[2]
			local jname2 : label (`by') `junk' 16
			local name1 = udsubstr(`"`jname1'"',1,16)
			local name2 = udsubstr(`"`jname2'"',1,16)
		}
		di _n in gr "Two-sample Kolmogorov–Smirnov test " /* 
			*/ "for equality of distribution functions" _n 
	}
	else {
		/* exp contains evaluation of expression */
		if "`exact'"!="" { error 198 } 
		capture assert (`exp')< . if `touse'
		if _rc {
			di in red "=exp evaluate to missing in " /*
			*/ "some or all observations"
			di in red "=exp should evaluate to cumulative " /*
			*/ "density evaluated at `varlist'"
			exit 450
		}
		local colexp = max(1, ceil((39-strlen("`o_exp'"))/2)+1)
		tempvar vn 
		qui gen double `vn' = `varlist' if `touse'
		local name1 = abbrev("`varlist'", 12)
		local name2 "Cumulative"
		di _n in gr "One-sample Kolmogorov–Smirnov test " /*
			*/ "against theoretical distribution" _n  /*
			*/ _col(`colexp') "`o_exp'" _n
	}

	tempvar G1 G2 V
	tempname m k Dm Dp D
	quietly {
		gen double `G1' = (`vn'< .) if `touse'
		replace `G1' = `weights'*(`G1'!=0) // NEW - use weights instead of 1
		gen double `V' = `exp' if `touse'
		if "`by'"!="" {
			gen double `G2' = `weights' if `V'< . // NEW - use weights instead of 1
			replace `V' = `vn' if `V'>=.
			sort `V'
			replace `G1' = sum(`G1')
			scalar `m' = `G1'[_N]
			replace `G1' = `G1' / `m'
			replace `G2' = sum(`G2')
			scalar `k' = `G2'[_N]
			replace `G2' = `G2' / `k'
			replace `G1' = `G1' - `G2'
			by `V': replace `G1' = `G1'[_N]
		}
		else {
			sort `vn'
			replace `G1' = sum(`G1')
			scalar `m' = `G1'[_N]
			replace `G1' = `G1' / `m' - `V'
			scalar `k' = 10e6
		}
		summarize `G1', meanonly
		scalar `Dm' = r(min)
		scalar `Dp' = r(max)
	
		if "`by'"=="" { scalar `Dm' = `Dm' - 1/`m' }
		scalar `D' = max(-`Dm',`Dp')
		
		// calculate p-value 
		tempname A Pm Pp Z2 P Pc Pz
		
		scalar `A' = `n1'*`n2'/(`n1'+`n2') // NEW - include effective sample sizes
		
		scalar `Pm' = exp(-2*`A'*((`Dm')^2))
		scalar `Pp' = exp(-2*`A'*((`Dp')^2))
		
		///// NEW - Use the asymptotic Q-function to approximate the 2-sided P-value /////
		scalar lambda = max((sqrt(`A')+0.12+.11/sqrt(`A'))*`D',0)
		local sum = 0
		forvalues j = 1/101 {
			local sum = `sum' + 2*(-1)^(`j'-1)*exp(-2*lambda^2*`j'^2)
		}
		scalar `P' = min(max(`sum',0), 1)
		scalar `Pc' = `P'
		
		tempname den eps
		if `P'<1 {
			if ("`exact'"!="") {
				local mk = `m'+`k'
				scalar `den' = /*
				*/ lngamma(`mk'+1)-lngamma(`m'+1)-lngamma(`k'+1)
				local i 0
				replace `G1' = 1 if _n==1
				scalar `eps' = 1/`m'/`k'/2
				while (`i'<`mk') {
					replace `G2' = `G1'[_n-1]+`G1'[_n]
					local i = `i' + 1
					replace `G2' = 1 if _n==1 | _n==`i'+1
					replace `G2' = 0 if abs((_n-1)/`m' ///
				-(`i'-_n+1)/`k')>=`D'-`eps' | _n>`i'+1
					replace `G1' = `G2'
				}
				if (`G1'[`m'+1] == 0) {
					scalar `Pc' = 1
				}
				else {	
					scalar `Pc' = log(`G1'[`m'+1]) - `den'
					scalar `Pc' = -expm1(`Pc')
				}	
			}
			else if `P' ~= 0 {
				scalar `Pz' = invnorm(`P')+1.04/min(`m',`k') /*
					*/ +2.09/max(`m',`k')-1.35/sqrt(`A')
				scalar `Pc' = normprob(`Pz')
			}
		}
	}
	if "`exact'"=="" {
		if (_caller()<14 | "`corrected'" != "") {
			local lbl "Corrected" 
		}
	}	
	else {
		local lbl "    Exact"
	}	
	di in gr "Smaller group             D     p-value  `lbl'"
	if ( "`exact'"!="" | (_caller()<14 | "`corrected'" != "")) {
		di in smcl in gr "{hline 50}"
	}
	else {
		di in smcl in gr "{hline 39}"
	}	
	di in gr "`name1'" _col(20) in ye %8.4f `Dp' "   " %9.3f `Pp'
	di in gr "`name2'" _col(20) in ye %8.4f `Dm' "   " %9.3f `Pm'
	if ( "`exact'"!="" | (_caller()<14 | "`corrected'" != "")) {
		di in gr "Combined K-S" _col(20) in ye %8.4f ///
			`D' "   " %9.3f `P' %11.3f `Pc'
	}
	else {
		di in gr "Combined K-S" _col(20) in ye %8.4f ///
			`D' "   " %9.3f `P' 
	}		
	/* double save in S_# and r() */
	ret local group1 "`name1'"
	ret local group2 "`name2'"
	ret scalar D_1 = `Dp'
	ret scalar D_2 = `Dm'
	ret scalar p_1 = `Pp'
	ret scalar p_2 = `Pm'
	ret scalar D   = `D'
	ret scalar p   = `P'
	if "`exact'"=="" { 
		if (_caller()<14 | "`corrected'" != "") {
			ret scalar p_cor   = `Pc' 
		}
		else {
			ret hidden scalar p_cor = `Pc'
		}	
	}	
	else    {    
		ret scalar p_exact = `Pc' 
	}	
	global S_1 "`name1'"
	global S_2 = `Dp'
	global S_3 = `Pp'
	global S_4 "`name2'"
	global S_5 = `Dm'
	global S_6 = `Pm'
	global S_8 = `D'
	global S_9 = `P'
	global S_10 = `Pc'

        qui count if `touse'
        local total = r(N)

        tempvar tag
        qui bys `touse' `varlist': gen byte `tag'= 1 if _n==1 & `touse'

        qui count if `tag'==1 & !missing(`v') & `touse'
        local unique = r(N)

        if `total' != `unique'{
		di
                if "`by'" == "" {
                        di as txt "Note: Ties exist in dataset; " _n ///
                        _col(7) "there are " string(`unique', "%9.0g") ///
			" unique values out of " string(`total', "%9.0g") ///
			" observations."
                }
                else {
                        di as txt "Note: Ties exist in combined " ///
			"dataset;" _n ///
                        _col(7) "there are " string(`unique', "%9.0g") ///
			" unique values out of " /// 
                        string(`total', "%9.0g") " observations."
                }
        }
end
