program define makewinnerb, rclass
// This function computes and draws subgroup winner/loser probability bounds conditional on subgroups of untreated outcomes.
// Implementation for bounds and confidence intervals is achieved using Chernozhukov, Lee, and Rosen (2013). 
// It also admits an instrument variable, results of which are conditioning on compliers.

// REQUIREMENT: need to install -moremata- by: ssc install moremata, replace

// When implementing this function, the user can increase nboot or decrease sparsity for higher precision at the cost of computation time. 

// By Jianfei Cao & Tetsuya Kaji
// Last updated: Nov 6, 2023

// structure of the program
// 1. stata main routines
// 2. stata subroutines
// 3. mata subroutines

////////// STATA MIAN ROUTINES //////////

syntax varlist [if] [in] [aw fw pw iw/] ///
	[,	LEFTtail /// draw lower quantile subgroups
		RIGHTtail /// draw upper quantile subgroups
		nboot(string) /// number of bootstrap replication in calculate confidence interval
		sparsity(string) /// set the gap width in the grid of the output graph
		instrument(varname) /// set intrument variable in monotone compliance 
		graphoff /// suppress graph output
		cboff /// suppress calculation of confidence interval (to save computation time)
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
qui gen `_y' = `1' // outcome variable 
qui gen `_d' = `2' // treatment status 
if `"`exp'"' != "" {
	gen `_w' = `exp' // weights
}
else {
	gen `_w' = 1 // set weights to 1 if not specified
}
qui putmata `_y' `_d' `_w' if `touse', replace

// select direction; default is to produce left tail bound
if "`lefttail'" != "" & "`righttail'" != "" {
	di as err "options {bf:left} and {bf:right} are mutually exclusive"
	exit 198
}

// default sparsity is 0.01
if ("`sparsity'" == "") {
	mata: sparsity = .01
}
else {
	mata: sparsity = strtoreal(st_local("sparsity"))
}
mata: u = evalPoints((`_y',`_d',`_w'),sparsity) // evaluation points
if "`righttail'" != "" {
	// grid for upper quantile subgroups
	mata: interMat = u,J(rows(u),1,1)
}
else {
	// grid for lower quantile subgroups
	mata: interMat = J(rows(u),1,0),u
}
mata: intervalMat = select(interMat,interMat[.,1]:!=interMat[.,2])
mata: u = select(u,interMat[.,1]:!=interMat[.,2])

// translate Stata variable into Mata matrix
if ("`instrument'" != "") {
	// monotone compliance case
	mata: ifInstrument = 1
	tempvar _z
	gen `_z' = `instrument'
	qui putmata `_z' if `touse', replace
	mata: data = (`_y',`_d',`_w',`_z')
}
else {
	// other cases
	mata: ifInstrument = 0
	mata: data = (`_y',`_d',`_w')
}

// set number of bootrap replication in implementing CLR bounds
if ("`nboot'" == "") {
	mata: nBoot = 1000 // defaulty is 1000
}
else {
	mata: nBoot = strtoreal(st_local("nboot"))
}

// main computation
if ("`cboff'" == "") {
	mata: results = clrBound(data, intervalMat, ifInstrument, 1, nBoot,sparsity)'
	mata: bndEstimate = results[.,1..4]
	mata: bndCI = results[.,5..8]
}
else {
	mata: bndEstimate = clrBound(data, intervalMat, ifInstrument, 0, nBoot,sparsity)'
}
	
// draw bounds
if ("`graphoff'" == "") {

	if ("`cboff'" == "") {
		// when confidence intervals are calculated 
		
		// temporarily delete dataset
		tempfile df
		qui save `df', replace
		clear 
		
		// translate Mata matrices back to Stata variables
		getmata u = u, replace
		getmata (winL winU losL losU) = bndEstimate, replace
		getmata (winLCI winUCI losLCI losUCI) = bndCI, replace
		qui replace losL = 1-losL
		qui replace losLCI = 1-losLCI
		gen zero = 0
		gen one = 1

		// draw
		graph twoway line winL u, ysc(r(0,1)) || /// 
			line losL u, ysc(r(0,1)) || ///
			rarea zero winLCI u, color(navy*.2) lwidth(none)|| ///
			(rarea losLCI one u, color(maroon*.2) lwidth(none)), ///
			legend(label(1 "winner lower bound") ///
				label(2 "1-(loser lower bound)") label(3 "CI lower bound") label(4 "CI upper bound")) ///
			xtitle("Rank") ///

		use `df', clear
	}
	else {
		// when confidence intervals are suppressed
		
		// temporarily delete dataset
		tempfile df
		qui save `df', replace
		clear 
		
		// translate Mata matrices back to Stata variables
		getmata u = u, replace
		getmata (winL winU losL losU) = bndEstimate, replace
		qui replace losL = 1-losL
		gen zero = 0
		gen one = 1

		// draw
		graph twoway line winL u, ysc(r(0,1)) || /// 
			line losL u, ysc(r(0,1)) ///
			legend(label(1 "winner lower bound") ///
			label(2 "1-(loser lower bound)")) ///
			xtitle("Rank") ///

		use `df', clear // get dataset back 
	}
}

// output data to result.csv
if ("`out2csv'" != "") {
	// temporarily delete dataset
	tempfile df
	qui save `df', replace
	clear
	
	// translate Mata matrices back to Stata variables 
	getmata (a b) = intervalMat, replace 
	getmata (winL winU losL losU) = bndEstimate, replace
	getmata (winLCI winUCI losLCI losUCI) = bndCI, replace
		
	capture mkdir "results"
	outsheet a b winL winU losL losU winLCI winUCI losLCI losUCI using results/`out2csv'.csv, comma replace
		
	use `df', clear // get dataset back 
}

// stored results
mata: st_local("l", strofreal(rows(u)))
mata: st_local("nBoot", strofreal(nBoot))
mata: st_local("sparsity", strofreal(sparsity))
return scalar l = `l' // number of evaluation points on the grid
return scalar sparsity = `sparsity' // gap width of the grid 
return scalar nBoot = `nBoot' // number of bootstrap replications
mata{
	st_matrix("a", intervalMat[.,1])
	st_matrix("b", intervalMat[.,2])
	st_matrix("winL", bndEstimate[.,1])
	st_matrix("winU", bndEstimate[.,2])
	st_matrix("losL", bndEstimate[.,3])
	st_matrix("losU", bndEstimate[.,4])
}
return matrix a a // all left ends of subgroups
return matrix b b // all right ends of subgroups
return matrix winL winL // lower winner bounds
return matrix winU winU // upper winner bounds 
return matrix losL losL // lower loser bounds
return matrix losU losU // upper loser bounds 

if ("`cboff'" == ""){
	mata: st_matrix("winLCI", bndCI[.,1])
	mata: st_matrix("winUCI", bndCI[.,2])
	mata: st_matrix("losLCI", bndCI[.,3])
	mata: st_matrix("losUCI", bndCI[.,4])
	return matrix winLCI winLCI // lower bounds of winner confidence intervals
	return matrix winUCI winUCI // upper bounds of winner confidence intervals
	return matrix losLCI losLCI // lower bounds of loser confidence intervals
	return matrix losUCI losUCI // upper bounds of loser confidence intervals
}

end


////////// MATA SUBROUTINES //////////

mata:

// calculates winner/loser bounds using Chernozhukov, Lee, and Rosen (2013)
function clrBound(data,interMat,ifInstrument,ifCI,nBoot,sparsity) {
			
	// data cleaning
	Y = data[.,1]
	D = data[.,2]
	W = data[.,3]

	// calculate empirical cdf
	if (ifInstrument==1) {
		// monotone compliance case
		Z = data[.,4]
		ecdf1 = ecdf0 = .
		makeecdfMC(ecdf1,ecdf0,data)
	}
	else {
		// other cases
		Y1 = select(Y,D:==1)	
		W1 = select(W,D:==1)
		ecdf1 = makeecdf(Y1,W1)
		Y0 = select(Y,D:==0)
		W0 = select(W,D:==0)
		ecdf0 = makeecdf(Y0,W0)
	}
	y1 = ecdf1[.,2]
	u1 = ecdf1[.,1]
	y0 = ecdf0[.,2]
	u0 = ecdf0[.,1]
	
	// evaluation points to add 
	x = uniqrows((interMat[.,1]',interMat[.,2]')')
	Q0 = evalEcdf(u0,y0,x,"Left")
	temp = sparseU(u0,sparsity)
	x = x\temp[.,1]
	Q0 = Q0\y0[temp[.,2],1]

	// main loop - i=1 for winner bounds, i=2 for loser bounds
	for (i=1; i<=2; i++) {
		
		// calculate true F1Q0
		if (i==1)	F1Q0 = evalEcdf(y1,u1,Q0)
		else 		F1Q0 = evalEcdf(y1,u1,Q0,"Left")
				
		// bootstrap F1Q0
		if (i==1)	bF1Q0 = f1q0Bootstrap(data,x,ifInstrument,nBoot)
		else 		bF1Q0 = f1q0Bootstrap(data,x,ifInstrument,nBoot,"Left")
		
		se = sqrt(diagonal(variance(bF1Q0'))) // bootstrap standard error 
		C = editmissing(correlation(bF1Q0'),0) // bootstrap coorelation matrix 
		
		// adaptive inequality selection
		n1 = rows(Y1)
		n0 = rows(Y0)
		gamma = 1-.1/log(rowmin((n1,n0))/2)
		a = interMat[.,1]'
		b = interMat[.,2]'
		V = ((a#J(rows(x),1,1)):<=(x#J(1,cols(a),1))):*((x#J(1,cols(a),1)):<=(b#J(rows(x),1,1)))
		k = normmaxq(C,V,gamma)
		
		// Preliminary set estimator
		lhs1 = ((x:-F1Q0)#J(1,cols(a),1))-(a#J(rows(x),1,1))
		nV1 = J(rows(V),cols(V),-100*max(abs(lhs1-se*k)))
		rhs1 = colmax(nV1:*(1:-V)+lhs1-se*k)#J(rows(V),1,1)-2*se*k
		V1 = lhs1:>=rhs1
		lhs2 = ((F1Q0-x:-1)#J(1,cols(b),1)):+(b#J(rows(x),1,1))
		nV2 = J(rows(V),cols(V),-100*max(abs(lhs2-se*k)))
		rhs2 = colmax(nV2:*(1:-V)+lhs2-se*k)#J(rows(V),1,1)-2*se*k
		V2 = lhs2:>=rhs2
		
		// Principal estimator
		k1 = normmaxq(C,V:*V1,.5)
		k2 = normmaxq(C,V:*V2,.5)
		sup1 = colmax(V:*(lhs1-se*k1):*(lhs1-se*k1:>=0))
		sup2 = colmax(V:*(lhs2-se*k2):*(lhs2-se*k2:>=0))
		
		// Convert to conditional probability bounds
		if (i==1) {
			winL = colmax(colmin(sup1:/(b-a)\J(1,cols(a),1))\J(1,cols(a),0))
			winU = colmax(colmin(1:-sup2:/(b-a)\J(1,cols(a),1))\J(1,cols(a),0))
		}
		else {
			losL = colmax(colmin(sup2:/(b-a)\J(1,cols(a),1))\J(1,cols(a),0))
			losU = colmax(colmin(1:-sup1:/(b-a)\J(1,cols(a),1))\J(1,cols(a),0))
		}

		if (ifCI==1) {
			// Principal estimator
			k1 = normmaxq(C,V:*V1,.95)
			k2 = normmaxq(C,V:*V2,.95)
			sup1 = colmax(V:*(lhs1-se*k1):*(lhs1-se*k1:>=0))
			sup2 = colmax(V:*(lhs2-se*k2):*(lhs2-se*k2:>=0))
			
			// Convert to conditional probability bounds
			if (i==1) {
				winLCI = colmax(colmin(sup1:/(b-a)\J(1,cols(a),1))\J(1,cols(a),0))
				winUCI = colmax(colmin(1:-sup2:/(b-a)\J(1,cols(a),1))\J(1,cols(a),0))
			}
			else {
				losLCI = colmax(colmin(sup2:/(b-a)\J(1,cols(a),1))\J(1,cols(a),0))
				losUCI = colmax(colmin(1:-sup1:/(b-a)\J(1,cols(a),1))\J(1,cols(a),0))
			}
		}
		
	}
	
	if (ifCI==1) {
		results = winL\winU\losL\losU\winLCI\winUCI\losLCI\losUCI
	}
	else {
		results = winL\winU\losL\losU
	}
	
	return(results)
	
}

// simulates F1(Q0()) function on evalations points 
function f1q0Bootstrap(data,x,ifInstrument,nBoot,|ifLeft) {
	
	// data cleaning 
	Y = data[.,1]
	D = data[.,2]
	W = data[.,3]:/sum(data[.,3]) // normalized weights
	n = rows(Y)
	Y1 = select(Y,D:==1)
	Y0 = select(Y,D:==0)
	bF1Q0 = J(rows(x),nBoot,0) // pre-alocate
		
	// main loop 
	for (i=1; i<=nBoot; i++) {
				
		// Poissonization - weights are drawn from Poisson(n*weight)
		wPoisson = rpoisson(1,1,n*W)

		if (ifInstrument==1) {
			// monotone compliance case
			ecdf1Boot = ecdf0Boot = .
			dataBoot = (Y,D,wPoisson,data[.,4])
			makeecdfMC(ecdf1Boot,ecdf0Boot,dataBoot)
		}
		else {
			// other cases 
			ecdf1Boot = makeecdf(Y1,select(wPoisson,D:==1))
			ecdf0Boot = makeecdf(Y0,select(wPoisson,D:==0))
		}
		y1 = ecdf1Boot[.,2]
		u1 = ecdf1Boot[.,1]
		y0 = ecdf0Boot[.,2]
		u0 = ecdf0Boot[.,1]

		Q0 = evalEcdf(u0,y0,x,"Left")
		if (args()==4)	bF1Q0[.,i] = evalEcdf(y1,u1,Q0)
		else			bF1Q0[.,i] = evalEcdf(y1,u1,Q0,"Left")
		
	}
		
	return(bF1Q0)
	
}

// calculates evaluation points needed for the grid 
function evalPoints(data,sparsity) {
	
	// data cleaning 
	Y = data[.,1]
	D = data[.,2]
	W = data[.,3]
	
	// calculate Y1 empircal cdf
	Y1 = select(Y,D:==1)
	W1 = select(W,D:==1)
	ecdf1 = makeecdf(Y1,W1)
	u1 = ecdf1[.,1]

	// calculate Y0 empircal cdf
	Y0 = select(Y,D:==0)
	W0 = select(W,D:==0)
	ecdf0 = makeecdf(Y0,W0)
	u0 = ecdf0[.,1]
	
	precision = .0000001
	u = sort(uniqrows((u0',(1:-u1)')')\precision\(1-precision),1)
	temp = sparseU(u,sparsity)
	u = temp[.,1]
	
	return(u)
}

// simulates joint normal distribution and calculates quantile of maximum 
function normmaxq(C,V,gamma) {

	// simulate normal with correlation C
	L = cholesky(C + (0.0000001)*I(rows(C)))
	simNorm = rnormal(500,rows(L),0,1)*L';

	// compute quantile of maximum of normals
	q = J(1,cols(V),0)
	for (i=1; i<=cols(q); i++) {
		q[1,i] = mm_quantile(rowmax(select(simNorm,V[.,i]')),1,gamma)
	}

	return(q)
	
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


