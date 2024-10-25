program define drawquantile
// This function draws and compares quantile functions of treated and control group outcomes.
// Weights are allowed. 
// It also admits an instrument variable, which draws quantiles conditional on compliers. 
// REQUIREMENT: need to install -moremata- by: ssc install moremata, replace

// By Jianfei Cao & Tetsuya Kaji
// Last updated: Oct 29, 2023

// structure of the program
// 1. stata main routines
// 2. mata subroutines


////////// STATA MAIN ROUTINES //////////

// initialization
syntax varlist [if] [in] [aw fw pw iw/] ///
	[, 	ymax(string) /// maximum value of y-axis
		ymin(string) /// minimum value of y-axis
		instrument(varname) /// instrument variable
		graphoff /// turn off graph to speed up computation
	]
marksample touse // mark sample selected by if/in command
markout `touse' `by', strok
confirm numeric variable `varlist' // string not allowed

// data cleaning
tempvar _y _d _w // outcome treatment weights 
gen `_y' = `1' // 1st input as outcome
gen `_d' = `2' // 2nd input as treatment
if `"`exp'"' != "" {
	gen `_w' = `exp' // weights if specified
}
else {
	gen `_w' = 1 // weights = 1 if not specified
}
qui putmata `_y' `_d' `_w' if `touse', replace // put data into mata

// calculate empirical cdf
if ("`instrument'" != "") {
	// when instrument is specified, calculate ecdf conditional on compliers
	tempvar _z
	gen `_z' = `instrument'
	qui putmata `_z' if `touse', replace // put data into mata
	mata: ecdf1 = ecdf0 = . // initialize ecdf
	mata: makeecdfMC(ecdf1,ecdf0,(`_y',`_d',`_w',`_z'))
}
else {
	// when instrument is not specified, calculate unconditional ecdf 
	mata: ecdf1 = makeecdf(select(`_y',`_d':==1),select(`_w',`_d':==1))
	mata: ecdf0 = makeecdf(select(`_y',`_d':==0),select(`_w',`_d':==0))
}

// determine y limit 
if ("`ymax'" == "") {
	// if ymax is not specified, set ymax to the maximum of ecdf1 and ecdf0
	mata: ymaxMata = max((max(ecdf1[.,2]),max(ecdf0[.,2])))
	mata: st_numscalar("ymax",ymaxMata) // store ymax in stata
	global ymax = ymax
}
else {
	// if ymax is specified
	global ymax = real("`ymax'")
}
if ("`ymin'" == "") {
	// if ymin is not specified, set ymin to the minimum of ecdf1 and ecdf0
	mata: yminMata = min((min(ecdf1[.,2]),min(ecdf0[.,2])))
	mata: st_numscalar("ymin",yminMata) // store ymin in stata
	global ymin = ymin
}
else {
	// if ymin is specified
	global ymin = real("`ymin'")
}

if ("`graphoff'" == "") {
	// draw the plot when not suppressing graph
	tempfile df
	qui save `df', replace
	clear // temporiarily clear data
	getmata (u1 y1) = ecdf1, replace // get ecdf1 from mata
	getmata (u0 y0) = ecdf0, replace force // get ecdf0 from mata
	graph set window fontface "Times New Roman"
	twoway (line y1 u1 if y1<=$ymax & y1>=$ymin, lpattern(dash) legend(label(1 "Q1")) lcolor(teal)) || ///
		(line y0 u0 if y0>=$ymin & y0<=$ymax, legend(label(2 "Q0")) lcolor(gs6)), ///
		yscale(range($ymin $ymax))
	use `df', clear // restore data
}

end


////////// MATA SUBROUTINES //////////
mata:

// make empirical cdf
function makeecdf(y,w) {
	
	// sort y and w
	dataSorted = sort((y,w),1)
	y = dataSorted[.,1]
	w = dataSorted[.,2]

	// obtain ecdf
	wCum = runningsum(w) // partial sum of weights
	yUniq = uniqrows(y) 
	ind = binarySearch(yUniq,y,"Last") // find the last index of each unique value
	uVec = wCum[ind,1] 
	uVec = uVec:/uVec[rows(uVec),1] // normalize
	
	// make ecdf drawable by specifying each turning point of step function
	uVec = uVec#(1\1) 
	uVec = 0\uVec[1..rows(uVec)-1,.]
	yVec = yUniq#(1\1)

	return((uVec,yVec))
	
}

// make empirical cdf of both treated and untreated compliers
function makeecdfMC(ecdf1,ecdf0,data) {
	
	// data cleaning
	y = data[.,1]
	d = data[.,2]
	w = data[.,3]
	z = data[.,4]

	in = z:*(1:-d) // index for never-takers
	ia = (1:-z):*d // index for always-takers
	
	ln = sum(in)/sum(z) // proportion of never-takers
	la = sum(ia)/sum(1:-z) // proportion of always-takers
	lc = 1-la-ln // proportion of compliers
	if (lc <= 0) _error("There are no compliers.")
	
	ecdfN = makeecdf(select(y,in:==1),select(w,in:==1)) // Never-taker cdf
	ecdfA = makeecdf(select(y,ia:==1),select(w,ia:==1)) // Always-taker cdf
	ecdfZ1 = makeecdf(select(y,z:==1),select(w,z:==1)) // Z=1 cdf
	ecdfZ0 = makeecdf(select(y,z:==0),select(w,z:==0)) // Z=0 cdf
	
	// if everyone is complier
	if (ln == 0 & la == 0) {
		ecdf1 = ecdfZ1
		ecdf0 = ecdfZ0
		return
	}
	
	// merge cdf of always-takers and never-takers
	Yna = ecdfN[.,2]\ecdfA[.,2]
	Yna = uniqrows(Yna)
	
	if (la > 0) {
		// If there are always-takers
		Una = la:*evalEcdf(ecdfA[.,2],ecdfA[.,1],Yna)
	}
	else {
		Una = J(rows(Yna),1,0)
	}
	
	if (ln > 0) {
		// If there are never-takers
		Una = Una+ln:*evalEcdf(ecdfN[.,2],ecdfN[.,1],Yna)
	}
	
	// make the low end 0 for eval_ecdf
	Yna = colmin(Yna)\Yna
	Una = 0\Una
	
	ecdf0 = constructMC(Yna,Una,lc,ecdfZ0)	// construct F0
	ecdf1 = constructMC(Yna,Una,lc,ecdfZ1)	// construct F1

}

// reusable function used to construct compliers ecdf
function constructMC(Yna,Una,lc,ecdfZ) {
	Yc = Yna\ecdfZ[.,2]
	Yc = uniqrows(Yc)
	UcU = evalEcdf(ecdfZ[.,2],ecdfZ[.,1],Yc)-evalEcdf(Yna,Una,Yc)
	UcU = UcU:/lc
	
	// rearrange to make it monotone in a mean-preserving way 
	dYc = Yc[2..rows(Yc),1]-Yc[1..rows(Yc)-1,1]
	temp = sort((UcU[1..rows(UcU)-1,1],dYc),(1,2))
	UcUs = temp[.,1]
	dYc = temp[.,2]
	UcUs = mm_cond(UcUs:<0,0,UcUs)
	UcUs = mm_cond(UcUs:>1,1,UcUs) // This is technically not mean-preserving, but should be negligible
	Yc = Yc[1,1]\(Yc[1,1]:+runningsum(dYc))

	// augment to make a plottable step function
	Uc = UcUs\1
	Uc = Uc#(1\1)
	Uc = 0\Uc[1..rows(Uc)-1,.]
	Yc = Yc#(1\1)
	
	return(uniqrows((Uc,Yc)))
}

// function that evaluate ecdf 
function evalEcdf(y,u,x,|ifLeft) {
	ind1 = J(rows(x),1,0)		
	if (args()==3) {
		ind = binarySearch(x,y,"Last")
	}
	else {
		ind = binarySearch(x,y)
	}
	return(u[ind,1])
}

// vectorized binary search
function binarySearch(A,B,|ifLast) {
// Syntax:
// L = binarySearch(A,B) returns an array containing the largest indices of elements of B that are smaller than A.
// L = binarySearch(A,B,"last") returns an array containing the smallest indices of elements of B that are larger than A.

	// initialization
	ub = rows(B)
	L = J(rows(A),1,1) // lower bounds
	U = J(rows(A),1,ub) // upper bounds

	// Handle A that exceeds max(B)
	L = mm_cond(A:>max(B),rows(B),L)
	
	// Logical index for remaining elements
	I = (B[1,.]:<=A):*(A:<=B[ub,1])
	
	// main loop
	while (sum(I)>0) {
		
		if (args()==2) {
			// ifLast is not specified
			M = floor((select(L,I)+select(U,I))/2)
			i1 = B[M,.]:<select(A,I)
		}
		else {
			// ifLast is specified
			M = ceil((select(L,I)+select(U,I))/2)
			i1 = B[M,.]:<=select(A,I)
		}
		i2 = I
		i2[selectindex(I),1] = i1
		if (sum(i2)>0) {
			L[selectindex(i2),1] = select(M,i1)
		}
		i2[selectindex(I),1] = 1:-i1
		if (sum(i2)>0) {
			U[selectindex(i2),1] = select(M,1:-i1)
		}
		I[selectindex(I),1] = ((select(L,I):+1) :< select(U,I))
		
	}

	// Adjust L if there are exact matches
	I = (B[U,.] :== A)
	L[selectindex(I),1] = U[selectindex(I),1]
	if (args()==2) L = mm_cond(A:==B[1,1],1,L)
	
	return(L)
}


end

