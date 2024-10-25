To use ste package locally, copy all files in ste-stata/ folder to the Stata personal folder. 
- Windows user: it is probably c:\ado\personal, but it might be someplace else.
- Mac user: it is probably ~/Documents/Stata/ado, but it might be someplace else. See Accessing personal ado-files to learn how Stata determines which folder to use.
- Unix user: it is ~/ado/personal (ado/personal in your home directory).
See https://www.stata.com/support/faqs/programming/personal-ado-directory/ for reference. 

In addition, this package requires -moremata- package. It can be installed by typing the following command in Stata:
ssc install moremata, replace

The package mainly contains 2 programs, makewfb and makewinnerb. 
- makewfb
- makewinnerb
*.ado file is the main program. *.sthlp contains help file. 

Examples/ folder contains replication files for all empirical results in the paper. Run example_tdj2015.do and example_bgh2006.do to obtain results. example_alo2009.do illustrates our methods in one-sided compliance cases to study treatment heterogeneity on compliers. 
