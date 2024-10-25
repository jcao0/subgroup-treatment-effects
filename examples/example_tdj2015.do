// version 18
// This script uses data from Tarozzi et al. (2015) to illustrate how to calculate bounds of subgroup treatment effect and corresponding confidence intervals. 

// Citation: Tarozzi, A., Desai, J., & Johnson, K. (2015). The impacts of microcredit: Evidence from Ethiopia. American Economic Journal: Applied Economics, 7(1), 54-89.


// initialization
// need to install -moremata- by:
// ssc install moremata, replace


*** TDJ2015 ***

clear all
set seed 7

timer clear 7
timer on 7


insheet using "data_tdj2015.csv", comma
keep if (time == "Endline") & (d_mf == 1 | d_non == 1)
capture mkdir "figures"

// Weighted Kolmogorov-Smirnov test
kstest2w cropnetsls, by(d_mf) weights(fwt_hh)
kstest2w anim_total, by(d_mf) weights(fwt_hh)

// Figure 1 - compare quantile functions
// (a)
drawquantile cropnetsls d_mf [aweight=fwt_hh], ymin(-1000) ymax(3000)
graph export figures/figure1a.pdf, replace

// (b)
drawquantile anim_total d_mf [aweight=fwt_hh], ymin(0) ymax(8000)
graph export figures/figure1b.pdf, replace


// Figure 2 -  STE bounds
// (a)
makewfb cropnetsls d_mf [aweight=fwt_hh], ymin(-1000) ymax(1000) nboot(10) out2csv(result)
graph export figures/figure2a.pdf, replace

// (b)
makewfb anim_total d_mf [aweight=fwt_hh], ymin(-1000) ymax(1000) nboot(10) 
graph export figures/figure2b.pdf, replace


// Figure 3 - welfare bounds
// (a) 
makewfb cropnetsls d_mf [aweight=fwt_hh], nrmloff ymin(0) nboot(10) 
graph export figures/figure3a.pdf, replace

// (b) 
makewfb anim_total d_mf [aweight=fwt_hh], nrmloff ymin(0) nboot(10)
graph export figures/figure3b.pdf, replace


// Figure 4 - non-utilitarian welfare
// (a) 
makewfb cropnetsls d_mf [aweight=fwt_hh], nrmloff ymin(0) custom nboot(10)
graph export figures/figure4a.pdf, replace

// (b)
makewfb anim_total d_mf [aweight=fwt_hh], nrmloff ymin(0) custom nboot(10)
graph export figures/figure4b.pdf, replace


timer off 7
timer list 7

