// version 18
// This script uses data from Angrist, Lang, and Oreopoulos (2009) to illustrate how to use IV to calculate subgroup treatment effect and subgroup winner probability for compliers. 

// Citation: Angrist, J., Lang, D., & Oreopoulos, P. (2009). Incentives and services for college achievement: Evidence from a randomized trial. American Economic Journal: Applied Economics, 1(1), 136-163.

// initialization
// need to install -moremata- by:
// ssc install moremata, replace


*** ALO2009 ***

clear all
set seed 7

timer clear 7
timer on 7

insheet using "data_alo2009.csv", comma
capture mkdir "figures"

drawquantile outcome treatment if gender == 0, instrument(eligibility)
graph export figures/oscQuantile.pdf, replace

makewfb outcome treatment if gender == 0, instrument(eligibility) nboot(10)
graph export figures/oscWfb.pdf, replace

makewinnerb outcome treatment if gender == 0, nboot(10) 
graph export figures/oscWinnerBound.pdf, replace

timer off 7
timer list 7

