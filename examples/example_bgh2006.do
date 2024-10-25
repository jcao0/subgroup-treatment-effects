// version 18
// This script uses data from Bloom et al. (1997) to illustrate how to calculate bounds of subgroup winner probability and corresponding confidence intervals. 

// Citation: Bloom, H. S., Orr, L. L., Bell, S. H., Cave, G., Doolittle, F., Lin, W., & Bos, J. M. (1997). The benefits and costs of JTPA Title II-A programs: Key findings from the National Job Training Partnership Act study. Journal of human resources, 549-576.

// initialization
// need to install -moremata- by:
// ssc install moremata, replace


*** BGH2006 ***

clear all
set seed 7


timer clear 7
timer on 7


insheet using "data_bgh2006.csv", comma
capture mkdir "figures"


// Figure 6
// (a)
drawquantile transq e if qtr == 4 [aweight = pscorewt], ymin(0) ymax(8000)
graph export figures/figure6a.pdf, replace

// (b)
makewinnerb transq e if qtr == 4 [aweight=pscorewt], nboot(10) out2csv(result)
graph export figures/figure6b.pdf, replace

// (c)
makewinnerb transq e if qtr == 4 [aweight=pscorewt], right nboot(10)
graph export figures/figure6c.pdf, replace

// (d)
drawquantile transq e if qtr == 12 [aweight = pscorewt], ymin(0) ymax(8000)
graph export figures/figure6d.pdf, replace

// (e)
makewinnerb transq e if qtr == 12 [aweight=pscorewt], nboot(10)
graph export figures/figure6e.pdf, replace

// (f)
makewinnerb transq e if qtr == 12 [aweight=pscorewt], righttail nboot(10)
graph export figures/figure6f.pdf, replace


// Figure 7
// (a)
drawquantile ernq e if qtr == 4 [aweight = pscorewt], ymin(0) ymax(8000)
graph export figures/figure7a.pdf, replace

// (b)
makewinnerb ernq e if qtr == 4 [aweight=pscorewt], nboot(10)
graph export figures/figure7b.pdf, replace

// (c)
makewinnerb ernq e if qtr == 4 [aweight=pscorewt], right nboot(10)
graph export figures/figure7c.pdf, replace

// (d)
drawquantile ernq e if qtr == 12 [aweight = pscorewt], ymin(0) ymax(8000) 
graph export figures/figure7d.pdf, replace

// (e)
makewinnerb ernq e if qtr == 12 [aweight=pscorewt], nboot(10)
graph export figures/figure7e.pdf, replace

// (f)
makewinnerb ernq e if qtr == 12 [aweight=pscorewt], righttail nboot(10)
graph export figures/figure7f.pdf, replace


timer off 7
timer list 7
