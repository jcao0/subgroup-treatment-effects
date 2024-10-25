program define makewfb, rclass
// This function computes and draws treatment effects/welfare change bounds conditional on subgroups of untreated outcomes.
// The main program integrates the difference (or transfermation of difference) between treatment and control quantile functions on appropriate intervals over the empirical measure. 
// It also admits an instrument variable, results of which are conditioning on compliers.

// REQUIREMENT: need to install -moremata- by: ssc install moremata, replace

// When implementing this function, the user can increase nboot or decrease sparsity for higher precision at the cost of computation time. 

// By Jianfei Cao & Tetsuya Kaji
// Last updated: Oct 29, 2023

// structure of the program
// 1. stata main routines
// 2. stata subroutines
// 3. mata subroutines

////////// STATA MIAN ROUTINES //////////

syntax varlist [if] [in] [aw fw pw iw/] ///
	[, 	te /// draw subgroup treatment effects 
		custom /// draw subgroup welfare effects using customized welfare function
		LEFTtail /// draw lower quantile subgroups
		RIGHTtail /// draw upper quantile subgroups
		nrmloff /// suppress normalization of conditional expectation 
		nboot(string) /// number of bootstrap replication in calculate confidence interval
		sparsity(string) /// set the gap width in the grid of the output graph
		instrument(varname) /// set intrument variable in monotone compliance 
		cboff /// suppress calculation of confidence interval (to save computation time)
		graphoff /// suppress graph output
		ymax(string) /// set max of y axis
		ymin(string) /// set min of y axis
		out2csv(string) /// output all results to a csv file
		seed(string) /// set seed for replication 
	] 

marksample touse // mark samples selected by if and in command
markout `touse' `by', strok
confirm numeric variable `varlist' // string not allowed

if ("`seed'" != "") {
	set seed `seed'
}

// data cleaning
tempvar _y _d _w
gen `_y' = `1' // outcome variable
gen `_d' = `2' // treatment status 
if `"`exp'"' != "" {
	gen `_w' = `exp' // weights
}
else {
	gen `_w' = 1 // set weights to 1 if not specified
}
qui putmata `_y' `_d' `_w' if `touse', replace

// select welfare function; default welfare function is identity
if ("`te'" != "") + ("`custom'" != "") > 1 {
	di as err "options {bf:te} and {bf:custom} are mutually exclusive"
	exit 198
}
local te = 0
mata: te = 0
if ("`custom'" != "") {
	mata: wfFcn = &wfFcnCustom() // customized welfare function
}
else {
	local te = 1 
	mata: te = 1
	mata: wfFcn = .
}

// select direction; default is to produce left tail bound
mata: righttail = 0
if ("`righttail'" != "") {
	mata: righttail = 1
}
if "`lefttail'" != "" & "`righttail'" != "" {
	di as err "options {bf:left} and {bf:right} are mutually exclusive"
	exit 198
}

// calculate empircal cdf of treated and control groups 
if ("`instrument'" != "") {
	// monotone compliance 
	mata: ifInstrument = 1
	tempvar _z // instrument variable
	gen `_z' = `instrument'
	qui putmata `_z' if `touse', replace
	mata: ecdf1 = ecdf0 = .
	mata: data = `_y',`_d',`_w',`_z'
	mata: makeecdfMC(ecdf1,ecdf0,data)
}
else {
	mata: ifInstrument = 0
	mata: data = `_y',`_d',`_w'
	mata: ecdf1 = makeecdf(select(`_y',`_d':==1),select(`_w',`_d':==1))
	mata: ecdf0 = makeecdf(select(`_y',`_d':==0),select(`_w',`_d':==0))
}

// default sparsity is .01; use sparsity=0 for finest possible grid
if ("`sparsity'" == "") {
	mata: sparsity = .01
}
else {
	mata: sparsity = strtoreal(st_local("sparsity"))
}

// compute the set of evaluation points on the grid
mata{
	precision = .0000001
	if (sparsity != 0) {
		// use a equal-gap grid according to the specifed sparsity 
		uVec = range(0,1,sparsity)
		uVec = mm_cond(uVec:==0,precision,uVec) // avoid 0
		if (righttail == 1) {
			intervalMat = uVec,J(rows(uVec),1,1) 
		}
		else {
			intervalMat = J(rows(uVec),1,0),uVec
		}
	}
	else {
		// sparsity=0 case, ie, finest possible grid 
		// use the discreteness of ecdf to calculate all evaluation points 
		// range(0,1,.1) is added to avoid to deal with constant Q1 Q0
		u1 = ecdf1[.,1]
		u0 = ecdf0[.,1]
		if (righttail == 1) {
			// upper quantile subgroups
			uL = sort(uniqrows(u0\(1:-u1)),1)
			uU = sort(uniqrows(u0\u1),1)
			uVec = sort(uniqrows(uU\uL\range(0,1,.1)),1)
			temp = mm_cond(uVec:==1,1-precision,uVec)
			uVec = temp[.,1]
			intervalMat = uVec,J(rows(uVec),1,1)
		}
		else {
			// lower quantile subgroups
			uL = u0#J(1,rows(u1),1)+u1'#J(rows(u0),1,1)
			uL = sort(uniqrows(colshape(mm_cond(uL:>1,1,uL),1)),1)
			uU = (1:+u0#J(1,rows(u1),1))-u1'#J(rows(u0),1,1)
			uU = mm_cond(uU:<0,0,uU)
			uU = sort(uniqrows(colshape(mm_cond(uU:>1,1,uU),1)),1)
			uVec = sort(uniqrows(uU\uL\range(0,1,.1)),1)
			temp = mm_cond(uVec:==0,precision,uVec)
			uVec = temp[.,1]
			intervalMat = J(rows(uVec),1,0),uVec
		}
	}
}

// MAIN PROGRAM

// compute welfare bounds
if (`te' == 1) {
	// identity welfare function
	mata: wfBound = makeWfb(ecdf1,ecdf0,intervalMat)
}
else {
	// customized welfare function
	mata: wfBound = makeWfb(ecdf1,ecdf0,intervalMat,&*wfFcn)
}

// compute bootstrap confidence intervals
if ("`cboff'" == "") {
	if ("`nboot'" == "") {
		mata: nBoot = 100 // default number of bootstrap replications is 100; change to 1000 for higher precision
	}
	else {
		mata: nBoot = strtoreal(st_local("nboot"))
	}
	if (`te' == 1) {
		// identity welfare function
		mata: confBand = imBootstrap(data,intervalMat,ifInstrument,nBoot)
	}
	else {
		// customized welfare function
		mata: confBand = imBootstrap(data,intervalMat,ifInstrument,nBoot,&*wfFcn)
	}
}

// de-normalize, eg, welfare analysis 
if ("`nrmloff'" != "") {
	// de-normailze the welfare bounds
	mata: wfBound = wfBound:*(uVec#J(1,cols(wfBound),1))
	if ("`cboff'" == "") {
		// de-normalize confidence intervals
		mata: confBand = confBand:*(uVec#J(1,cols(confBand),1))
	}
}

// draw the plot
if ("`graphoff'" == "") {
	// determine y limit 
	// if not specifed, use max or min of welfare bounds or confidence intervals
	if ("`ymax'" == "") {
		if ("`cboff'" == ""){
			mata: ymaxMata = max(confBand)
		}
		else {
			mata: ymaxMata = max(wfBound)
		}
		mata: st_numscalar("ymax",ymaxMata)
		global ymax = ymax
	}
	else {
		global ymax = real("`ymax'") // when ymax is specified
	}
	if ("`ymin'" == "") {
		if ("`cboff'" == ""){
			mata: yminMata = min(confBand)/1.1
		}
		else {
			mata: yminMata = min(wfBound)/1.1
		}
		mata: st_numscalar("ymin",yminMata)
		global ymin = ymin
	}
	else {
		global ymin = real("`ymin'") // when ymin is specified
	}

	// determine whether to draw confidence intervals
	if ("`cboff'" == "") {
		drawCB // draw CI
	}
	else {
		drawNoCB // don't draw CI
	}
}

// output data to result.csv
if ("`out2csv'" != "") {

	tempfile df
	qui save `df', replace
	clear // temporarily delete the original dataset
	getmata (a b) = intervalMat, replace
	getmata (ub lb) = wfBound, replace force
		
	capture mkdir "results"

	if ("`cboff'" == ""){
		// case with CI
		getmata (ubCI lbCI) = confBand, replace force
		outsheet a b ub lb ubCI lbCI using results/`out2csv'.csv, comma replace
	}
	else {
		// case without CI
		outsheet uVec ub lb results/`out2csv'.csv , comma replace	
	}
	use `df', clear
}

// stored results
mata: st_local("l", strofreal(rows(uVec)))
mata: st_local("sparsity", strofreal(sparsity))
return scalar l = `l' // number of evaluation points on the grid
return scalar sparsity = `sparsity' // gap width of the grid 

mata{
	st_matrix("a", intervalMat[.,1])
	st_matrix("b", intervalMat[.,2])
	st_matrix("ub", wfBound[.,1])
	st_matrix("lb", wfBound[.,2])
}
return matrix a a // all left ends of subgroups
return matrix b b // all right ends of subgroups
return matrix ub ub // upper bounds
return matrix lb lb // lower bounds

if ("`cboff'" == ""){
	mata: st_local("nBoot", strofreal(nBoot))
	return scalar nBoot = `nBoot' // number of bootstrap replications
	
	mata: st_matrix("ubCI", confBand[.,1])
	mata: st_matrix("lbCI", confBand[.,2])
	return matrix ubCI ubCI // confidence interval upper bounds
	return matrix lbCI lbCI // confidence interval lower bounds
}

end


////////// STATA SUBROUTINES //////////

program define drawCB
// this program draws the welfare bounds along with confidence intervals

// temporarily delete the original dataset 
tempfile df
qui save `df', replace
clear

// translate Mata-format bounds into Stata dataset
getmata (ub lb) = wfBound
getmata (ubCI lbCI) = confBand, replace force
getmata uVec = uVec, replace force

// draw 
qui replace ubCI = $ymax if ubCI>$ymax
qui replace lbCI = $ymin if lbCI<$ymin
graph set window fontface "Times New Roman"
graph twoway rarea ubCI lbCI uVec, color(gs14) || ///
	line ub uVec if ub<=$ymax & ub>=$ymin, lpattern(dash) || ///
	line lb uVec if lb>=$ymin & lb<=$ymax, ///
	graphregion(color(white)) ///
	ysc(r($ymin,$ymax)) ///
	legend(cols(1) order( ///
	2 "upper bound" ///
	3 "lower bound" ///
	1 "Confidence interval" ///
	)) ///
	xtitle("Rank")
use `df', clear

end


program define drawNoCB
// this program draws the welfare bounds without confidence intervals

// temporarily delete the original dataset 
tempfile df
qui save `df', replace
clear

// translate Mata-format bounds into Stata dataset
getmata (ub lb) = wfBound
getmata uVec = uVec, replace force

// draw
graph set window fontface "Times New Roman"
graph twoway line lb uVec if lb>=$ymin & lb<=$ymax || ///
	line ub uVec if ub<=$ymax & ub>=$ymin, lpattern(dash) ///
	graphregion(color(white)) ///
	ysc(r($ymin,$ymax)) ///
	legend(label (1 "lower bound") label (2 "upper bound")) ///
	xtitle("Rank")
use `df', clear

end


////////// MATA SUBROUTINES //////////

mata:

// customize welfare function
// default example is a welfare function that penalizes loss over gain at a rate of 1.1:1
function wfFcnCustom(x) {
	return((x:>=0):*x:+(x:<0):*(1.1:*x))
}

// calculate welfare bounds
function makeWfb(ecdf1,ecdf0,intervalMat,|func) {
	// pre-define output matrix 
	upperBound = J(rows(intervalMat),1,0)
	lowerBound = J(rows(intervalMat),1,0)
	
	// main loop
	for (i=1; i<=rows(intervalMat); i++) {
		a = intervalMat[i,1]
		b = intervalMat[i,2]
		if (args()==3) {
			// identity welfare function
			upperBound[i,.] = shiftBound(ecdf1,ecdf0,a,b)
			lowerBound[i,.] = flipBound(ecdf1,ecdf0,a,b)
		} 
		else {
			// customized welfare function
			upperBound[i,.] = shiftBound(ecdf1,ecdf0,a,b,&*func)
			lowerBound[i,.] = flipBound(ecdf1,ecdf0,a,b,&*func)
		}
	}		
	return((upperBound,lowerBound))
}

// function that integrates f(Q1(1-b+u)-Q0(u)) over a to b
function shiftBound(q1,q0,a,b,|func) {
	
	// data cleaning
	u1 = q1[.,1]
	y1 = q1[.,2]
	u0 = q0[.,1]
	y0 = q0[.,2]

	// set relevant range in treatment and control quantile functions
	x1 = mm_cond(u1:<1-b+a,1-b+a,u1)
	x0 = mm_cond(u0:<a,a,u0)
	x0 = mm_cond(x0:>b,b,x0)

	// main calculation
	if (args()==4)  {
		// when welfare is identity
		q1Int = sum((x1[2..rows(x1),.]-x1[1..rows(x1)-1,.]):*y1[1..rows(y1)-1,.])
		q0Int = sum((x0[2..rows(x0),.]-x0[1..rows(x0)-1,.]):*y0[1..rows(y0)-1,.])	
		return((q1Int-q0Int)/(b-a))
	} 
	else {
		// customized welfare function
		
		// select relevant evaluation points 
		u1Temp = select(u1,u1[.,1]:>=1-b+a):-(1-b) // y1 grid
		u0Temp = select(u0, (u0[.,1]:>=a) :* (u0[.,1]:<=b)) // y0 grid
		y1y1Grid = select(y1,u1[.,1]:>=1-b+a)
		y0y0Grid = select(y0,(u0[.,1]:>=a) :* (u0[.,1]:<=b))
		
		// treatment part integrand	
		y1y0Grid = evalEcdf(u1,y1,(1-b):+u0Temp,"Left")
		q1ab = sort((u1Temp,y1y1Grid)\(u0Temp,y1y0Grid),(1,2))

		// control part integrand
		y0y1Grid = evalEcdf(u0,y0,u1Temp,"Left")
		q0ab = sort((u0Temp,y0y0Grid)\(u1Temp,y0y1Grid),(1,2))

		// integral
		y = (*func)(q1ab[.,2]-q0ab[.,2])
		x = q1ab[.,1]

		result = sum((x[2..rows(x),.]-x[1..rows(x)-1,.]):*y[1..rows(y)-1,.])/(b-a)
		return(result)		
	}
}

// function that integrates f(Q1(b-u)-Q0(u)) over a to b
function flipBound(q1,q0,a,b,|func) {

	// data cleaning
	u1 = q1[.,1]
	y1 = q1[.,2]
	u0 = q0[.,1]
	y0 = q0[.,2]

	// set relevant range in treatment and control quantile functions
	x1 = mm_cond(u1:>b-a,b-a,u1)
	x0 = mm_cond(u0:<a,a,u0)
	x0 = mm_cond(x0:>b,b,x0)
	
	if (args()==4) {
		// when welfare is identity
		q1Int = sum((x1[2..rows(x1),.]-x1[1..rows(x1)-1,.]):*y1[1..rows(y1)-1,.])
		q0Int = sum((x0[2..rows(x0),.]-x0[1..rows(x0)-1,.]):*y0[1..rows(y0)-1,.])
		return((q1Int-q0Int)/(b-a))
	} 
	else {
		// customized welfare function

		// select relevant evaluation points 
		u1Temp = b:-select(u1,u1[.,1]:<=b-a) // y1 grid
		u0Temp = select(u0, (u0[.,1]:>=a) :* (u0[.,1]:<=b)) // y0 grid
		y1y1Grid = select(y1,u1[.,1]:<=b-a)
		y0y0Grid = select(y0,(u0[.,1]:>=a) :* (u0[.,1]:<=b))
		
		// treatment part integrand
		y1y0Grid = evalEcdf(u1,y1,b:-u0Temp,"Left")
		q1ab = sort((u1Temp,y1y1Grid)\(u0Temp,y1y0Grid),(1,-2))

		// control part integrand 
		y0y1Grid = evalEcdf(u0,y0,u1Temp,"Left")
		q0ab = sort((u0Temp,y0y0Grid)\(u1Temp,y0y1Grid),(1,2))

		// integral
		y = (*func)(q1ab[.,2]-q0ab[.,2])
		x = q1ab[.,1]
		
		result = sum((x[2..rows(x),.]-x[1..rows(x)-1,.]):*y[1..rows(y)-1,.])/(b-a)
		return(result)
	}
}

// IM2004 bootstrap confidence band
function imBootstrap(data,intervalMat,ifInstrument,nBoot,|func) {

	// data cleaning
	Y = data[.,1] // outcome
	D = data[.,2] // treatment 
	W = data[.,3]
	W = W:/colsum(W) // normalized weights
	if (ifInstrument == 1) {
		Z = data[.,4] // instrument in monotone compliance
	}
	n = rows(Y)

	// main loop 
	sumBoot = 0
	sum2Boot = 0
	for (i=1; i<=nBoot; i++) {

		// Poissonization - weights are drawn from Poisson(n*weight)
		wPoisson = rpoisson(1,1,n*W)

		// calculate bootstrap ecdf; only weights are random across replications
		if (ifInstrument == 1) {
			ecdf1Boot = ecdf0Boot = .
			dataBoot = (Y,D,wPoisson,Z) 
			makeecdfMC(ecdf1Boot,ecdf0Boot,dataBoot)
		}
		else {
			if (sum(select(wPoisson,D:==1)):==0) {
				// the rare case when all weights are zero
				ecdf1Boot = makeecdf(select(Y,D:==1),select(W,D:==1))
			}
			else {
				ecdf1Boot = makeecdf(select(Y,D:==1),select(wPoisson,D:==1))
			}
			if (sum(select(wPoisson,D:==0)):==0) {
				// the rare case when all weights are zero
				ecdf0Boot = makeecdf(select(Y,D:==0),select(W,D:==0))
			}
			else {
				ecdf0Boot = makeecdf(select(Y,D:==0),select(wPoisson,D:==0))
			}		
		}
	
		// calculate bootstrap bounds
		if (args()==4) {
			wfBound = makeWfb(ecdf1Boot,ecdf0Boot,intervalMat)
		}
		else {	
			wfBound = makeWfb(ecdf1Boot,ecdf0Boot,intervalMat,&*func)
		}

		// record sum and sum of squares 
		sumBoot = sumBoot:+wfBound
		sum2Boot = sum2Boot:+wfBound:*wfBound

	}

	// recover mean and variance 
	meanBoot = sumBoot/nBoot
	varBoot = (sum2Boot-sumBoot:*sumBoot/nBoot)/(nBoot-1)
	varBoot = mm_cond(varBoot:<0,0,varBoot) // avoid negative variance
	sdBoot = sqrt(varBoot)
	
	// use Imbens and Manski (2004) bounds 
	imBound = im2004(meanBoot,sdBoot,.05)
	imBound = imBound:*(sdBoot:!=0)+meanBoot:*(sdBoot:==0) // deal with sd==0
	return(imBound)
	
}

// IM2004 Conficence interval of Imbens and Manski (2004)
function im2004(meanMat,seMat,alpha){
// Makes a confidence interval for a set-identified parameter. 
// This is different from a confidence interval for the identified set. 
// Relies on the normal approximation of the lower and upper bound estimators.

	C = J(rows(meanMat),1,0)
	maxSE = rowmax(seMat)
	delta = meanMat[.,1]:-meanMat[.,2]
	shift = delta:/maxSE

	for (i_u=1; i_u<=rows(meanMat); i_u++) {
		if (maxSE[i_u,1] == 0) C[i_u,1] = 0
		else {
			findRoot = mm_root(c=., &imSolver(), 0, 10, 0, 1000, shift[i_u,1],alpha)
			C[i_u,1] = c
		}
	}
	
	imLowerBound = meanMat[.,2]:-seMat[.,2]:*C
	imUpperBound = meanMat[.,1]:+seMat[.,1]:*C	
	return((imUpperBound,imLowerBound))
}

// function used in calculating IM2004 correction
function imSolver(c,shift,alpha) {
	return(normal(c+shift)-normal(-c)-1+alpha)
}

// make empirical cdf
function makeecdf(y,w) {
	
	// sort y and w
	dataSorted = sort((y,w),1)
	y = dataSorted[.,1]
	w = dataSorted[.,2]

	// obtain ecdf
	wCum = runningsum(w)
	yUniq = uniqrows(y)
	ind = binarySearch(yUniq,y,"Last")
	uVec = wCum[ind,1]
	uVec = uVec:/uVec[rows(uVec),1]
	
	// make ecdf drawable by adding all turning points
	uVec = uVec#(1\1)
	uVec = 0\uVec[1..rows(uVec)-1,.]
	yVec = yUniq#(1\1)

	return((uVec,yVec))
	
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

// sparsify evaluation points
function sparseU(u,sparsity) {
	if (sparsity == 0) {
		return((u,(1..rows(u))'))
	}
	else {
		precision = .0000001
		sparseVec = range(0,1,sparsity)
		sparseVec = mm_cond(sparseVec:==0,precision,sparseVec)
		sparseVec = mm_cond(sparseVec:==1,1-precision,sparseVec)
		ind = binarySearch(sparseVec,u)
		return((u[ind,1],ind))
	}
}

end










