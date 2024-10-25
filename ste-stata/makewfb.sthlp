{smcl}
{* *! version 1.1.10  18mar2021}{...}

{title:Title}

{p2colset 5 22 15 2}{...}
{p2col :{cmd:makewfb} {hline 2}}Draw subgroup welfare effects bounds introduced in Kaji and Cao (2023) {p_end}
{p2colreset}{...}

{marker syntax}{...}
{title:Syntax}

{p 8 15 2} {cmd:makewfb}
outcome treatment
{ifin} {it:{weight}} {cmd:,} [{help makewfb##options:options}] {p_end}

{marker options}{...}
{title:Options}

{marker opt_summary}{...}
{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Welfare function}
{synopt :{opt te}}(default) calculate subgroup treatment effect bounds E[Y1-Y0|U<b], where b ranges from 0 to 1 and U is "rank" of Y0 as in Y0=Q0(U);
in this case, welfare function is identity;
condition becomes U>b if direction is set to righttail{p_end}
{synopt :{opt custom}}other customized welfare function that user can specify in makewfb.ado; default is a precoded welfare function that penalizes loss at the rate of 1.1 to 1;
specifically, it calculates E[(Y1-Y0)1{Y1>Y0}+1.1*(Y1-Y0)1{Y1<=Y0}|U<b];
condition becomes U>b if direction is set to righttail 
{p_end}

{syntab:Direction}
{synopt :{opt left:tail}}(default) subgroup for lower quantiles, i.e., subgroup effect is conditional on [0,b] for b ranging from 0 to 1{p_end}
{synopt :{opt right:tail}}subgroup for upper quantiles, i.e., subgroup effect is conditional on [a,1] for a ranging from 0 to 1{p_end}

{syntab:Configuration}
{synopt :{opt nrmloff}}suppress normalization by subgroup proportion, e.g., reporting E[(Y1-Y0)1{U<b}] instead of E[Y1-Y0|U<b]{p_end}
{synopt :{opt nBoot(integer)}}number of bootstrap replications when calculating confidence intervals introduced by Imbens and Manski (2004); default is set to 100; use nBoot(1000) for higher precision at the cost of computation time{p_end}
{synopt :{opt sparsity(number)}}gap width in the grid that specifies subgroups, e.g., possible values of b when calculating E[Y1-Y0|U<b]; 
use sparsity(0) for finest possible grid, which is not recommended for large datasets; 
default is set to .01; 
try large values, say .1, if execution time is too long{p_end}
{synopt :{opt instrument(varname)}}monotone compliance - this option specifies an instrument variable (e.g., randomized eligibility) and produces results for compliers{p_end}
{synopt :{opt seed(number)}}user can set seed for reproducibility{p_end}

{syntab:Output}
{synopt :{opt graphoff}}suppress graph output{p_end}
{synopt :{opt cboff}}suppress confidence intervals, which greatly saves computation time{p_end}
{synopt :{opth ymax(number)}}maximam y axis{p_end}
{synopt :{opth ymin(number)}}minimum y axis{p_end}
{synopt :{opt out2csv(string)}}use out2csv(<filename>) to generate <filename>.csv that contains numerical results with following columns:{p_end}
{synopt :} a - lower end of range of Y0 rank{p_end}
{synopt :} b - upper end of range of Y0 rank{p_end}
{synopt :} ub - upper bound for subgroup treatment effect or welfare change{p_end}
{synopt :} lb - lower bound for subgroup treatment effect or welfare change{p_end}
{synopt :} ubCI - pointwise upper confidence interval{p_end}
{synopt :} lbCI - pointwise lower confidence interval{p_end}
{synoptline}

{marker description}{...}
{title:Description}

{pstd}{cmd:makewfb} calculates and draws the upper and lower sharp bounds of the treatment effects of a sequence of subgroups. 
This is, it asks questions such as what is the average effect of the poverty-reducation policy for the poorest 20% population. 
This function is useful assessing the heterogeneity of treatment effects across different subgroups and helping with welfare analysis.

{pstd}The function takes as input the outcome variable, the treatment indicator, and optionally the weights.
By default, the function draws estimates of bounds for E[f(Y1-Y0)|0<U<b] against b, 
where b ranges from 0 to 1, 
f is some welfare function (identity as the default), 
and U is the rank of Y0 (i.e., U=F0(Y0), where F0 is cdf of Y0). 
The user can choose to output a table that displays the upper and lower bounds of the welfare effects for each subgroup, along with the confidence intervals. 

{pstd}Estimation - The estimates of bounds for E[f(Y1-Y0)|0<U<b] are calculated by intergrating over the empirical measure. 

{pstd}Inference - The confidence intervals are calculated using Imben and Manski (2004) and are pointwise.
The confidence intervals are at 95% confidence level. 
The level is with respect to the probability of CI covering the set-identified parameter, instead of CI covering the identified set. 

{pstd}The function can deal with the monotone compliance case where the treatment effects or welfare change are with respect to the compliers.
Inputting an instrument variable is required in this setting. 


{marker contact}{...}
{title:Authors}

{pstd}Jianfei Cao{break}
Department of Economics, Northeastern University{break}
Email: {browse "mailto:j.cao@northeastern.edu":j.cao@northeastern.edu}
{p_end}

{pstd}Tetsuya Kaji{break}
University of Chicago Booth School of Business{break}
Email: {browse "mailto:tkaji@chicagobooth.edu":tkaji@chicagobooth.edu}
{p_end}

{marker citation}{...}
{title:Citation}

{pstd}Kaji, T., & Cao, J. (2023). Assessing Heterogeneity of Treatment Effects. arXiv preprint arXiv:2306.15048.

{pstd}Imbens, G. W., & Manski, C. F. (2004). Confidence intervals for partially identified parameters. Econometrica, 72(6), 1845-1857.

{marker support}{...}
{title:Support}

{pstd}{cmd:makewfb} requires the {cmd:moremata} package.

{marker examples}{...}
{title:Examples}

{pstd}Load dataset from Cattaneo (2010) Journal of Econometrics 155: 138–154, which studies the effect of maternal smoking intensity during pregnancy on birth weight{p_end}
{phang2}{cmd:. webuse cattaneo2}{p_end}

{pstd}Compare quantile functions of infant birthweight between treated and control group, where the treatment variable is whether mother smokes; Q1 is for smoking and Q0 for not smoking{p_end}
{phang2}{cmd:. drawquantile bweight mbsmoke}{p_end}

{pstd}Draw upper and lower bounds of subgroup treamtment effects E[Y1-Y0|U<b] against b, along with the confidence intervals{p_end}
{phang2}{cmd:. makewfb bweight mbsmoke}{p_end}

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:makewfb} stores the following in {cmd:r()}:

{synoptset 15 tabbed}{...}
{p2col 5 15 19 2: Scalars}{p_end}
{synopt:{cmd:r(l)}}number of evaluation points in quantile grid{p_end}
{synopt:{cmd:r(nBoot)}}number of bootstrap replications when calculating confidence intervals{p_end}
{synopt:{cmd:r(sparsity)}}gap width in quantile grid{p_end}

{synoptset 15 tabbed}{...}
{p2col 5 15 19 2: Matrices}{p_end}
{synopt:{cmd:r(a)}}lower bounds of subgroup ranges{p_end}
{synopt:{cmd:r(b)}}upper bounds of subgroup ranges{p_end}
{synopt:{cmd:r(ub)}}upper bound of subgroup welfare effects{p_end}
{synopt:{cmd:r(lb)}}lower bound of subgroup welfare effects{p_end}
{synopt:{cmd:r(ubCI)}}upper bound of confidence intervals{p_end}
{synopt:{cmd:r(lbCI)}}lower bound of confidence intervals{p_end}
{p2colreset}{...}





