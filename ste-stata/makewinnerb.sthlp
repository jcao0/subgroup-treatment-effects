{smcl}
{* *! version 1.1.10  18mar2021}{...}

{title:Title}

{p2colset 5 22 15 2}{...}
{p2col :{cmd:makewinnerb} {hline 2}}Draw subgroup winner bounds introduced in Kaji and Cao (2023){p_end}
{p2colreset}{...}

{marker syntax}{...}
{title:Syntax}

{p 8 15 2} {cmd:makewinnerb}
outcome treatment
{ifin} {it:{weight}} {cmd:,} [{help makewinnerb##options:options}] {p_end}

{marker options}{...}
{title:Options}

{marker opt_summary}{...}
{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Direction}
{synopt :{opt left:tail}}(default) subgroup for lower quantiles, i.e., subgroup probability is conditional on untreated outcome Y0 being within the lowest quantile [0,b], for various b{p_end}
{synopt :{opt right:tail}}subgroup for upper quantiles, i.e., Y0 within the highest quantile [a,1], for various a{p_end}

{syntab:Configuration}
{synopt :{opt nBoot(integer)}}number of bootstrap replications when calculating bounds and confidence intervals; default is set to 1000{p_end}
{synopt :{opt sparsity(number)}}gap width in quantile grid; use sparsity(0) for finest possible grid, not recommended for large datasets; default is set to 0.01{p_end}
{synopt :{opt instrument(varname)}}monotone compliance - this option specifies an instrument variable (e.g., randomized eligibility) and produces results for compliers{p_end}
{synopt :{opt seed(number)}}user can set seed for reproducibility{p_end}

{syntab:Output}
{synopt :{opt graphoff}}suppress graph output{p_end}
{synopt :{opt cboff}}suppress confidence intervals, which greatly saves computation time{p_end}
{synopt :{opth ymax(number)}}maximam yaxis{p_end}
{synopt :{opth ymin(number)}}minimum yaxis{p_end}
{synopt :{opt out2csv(string)}}use out2csv({it:filename}) to generate {it:filename}.csv that contains results with following columns:{p_end}
{synopt :} a - lower end of range of Y0 rank{p_end}
{synopt :} b - upper end of range of Y0 rank{p_end}
{synopt :} winL - lower bounds of winner proportions{p_end}
{synopt :} winU - upper bounds of winner proportions{p_end}
{synopt :} losL - lower bounds of loser proportions{p_end}
{synopt :} losU - upper bounds of loser proportions{p_end}
{synopt :} winLCI - lower bounds of confidence intervals of winner proportions{p_end}
{synopt :} winUCI - upper bounds of confidence intervals of winner proportions{p_end}
{synopt :} losLCI - lower bounds of confidence intervals of loser proportions{p_end}
{synopt :} losUCI - upper bounds of confidence intervals of loser proportions{p_end}
{synoptline}

{marker description}{...}
{title:Description}

{pstd}{cmd:makewinnerb} This function calculates and draws the upper and lower bounds of the proportion of winners and losers for a treatment in a sequence of subgroups.
This is, it asks questions such as how many of the poorest 20% population really benefit from the povery-reduction policy. 
By shedding light on who benefit from the treatment, this function is useful assessing the heterogeneity of treatment effects across different subgroups and helping with welfare analysis.

{pstd}The function takes as input the outcome variable, the treatment indicator, and optionally the weights.
By default, the function draws 1) the lower bound for Pr(Y1>Y0|0<U<b) against b, 
where b ranges from 0 to 1 and U is the rank of Y0 (i.e., U=F_0(Y0)), 
2) the upper bound for 1-Pr(Y1<Y0|0<U<b) against b, which is close to the upper bound for Pr(Y1>Y0|0<U<b) except for the Y1=Y0 case. 
The user can choose to output a table that displays the upper and lower bounds of the winner and loser proportion for each subgroup, along with the confidence intervals. 

{pstd}Both the estimates and the confidence intervals are calculated using Chernozhukov, Lee, and Rosen (2013). 

{pstd}Estimation - The bound estimates are median unbiased, in the sense that probability of the set-identified parameter being smaller than the upper bound estimate (or larger than the lower bound estimate) is asymptotically 0.5. 

{pstd}Inference - The confidence intervals are at 95% confidence level, in the sense that probability of the set-identified parameter being smaller than the upper bound estimate (or larger than the lower bound estimate) is asymptotically 0.95. 

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

{pstd}Chernozhukov, V., Lee, S., & Rosen, A. M. (2013). Intersection bounds: Estimation and inference. Econometrica, 81(2), 667-737.

{pstd}Chernozhukov, V., Kim, W., Lee, S., & Rosen, A. M. (2015). Implementing intersection bounds in Stata. The Stata Journal, 15(1), 21-44.

{marker support}{...}
{title:Support}

{pstd}{cmd:makewinnerb} requires the {cmd:moremata} package.

{marker examples}{...}
{title:Examples}

{pstd}Load dataset from Cattaneo (2010) Journal of Econometrics 155: 138–154, which studies the effect of maternal smoking intensity during pregnancy on birth weight{p_end}
{phang2}{cmd:. webuse cattaneo2}{p_end}

{pstd}Compare quantile functions of infant birthweight between treated and control group, where the treatment variable is whether mother smokes; Q1 is for smoking and Q0 for not smoking{p_end}
{phang2}{cmd:. drawquantile bweight mbsmoke}{p_end}

{pstd}Draw upper and lower bounds of subgroup winner probability Pr(Y1>Y0|U<b) against b, along with the confidence intervals{p_end}
{phang2}{cmd:. makewinnerb bweight mbsmoke}{p_end}

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
{synopt:{cmd:r(winL)}}lower bounds of winner proportions{p_end}
{synopt:{cmd:r(winU)}}upper bounds of winner proportions{p_end}
{synopt:{cmd:r(losL)}}lower bounds of loser proportions{p_end}
{synopt:{cmd:r(losU)}}upper bounds of loser proportions{p_end}
{synopt:{cmd:r(winLCI)}}lower bounds of confidence intervals of winner proportions{p_end}
{synopt:{cmd:r(winUCI)}}upper bounds of confidence intervals of winner proportions{p_end}
{synopt:{cmd:r(losLCI)}}lower bounds of confidence intervals of loser proportions{p_end}
{synopt:{cmd:r(losUCI)}}upper bounds of confidence intervals of loser proportions{p_end}
{p2colreset}{...}





