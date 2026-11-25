/*
master.do runs the entire pipeline. sets up necessary globals, makes necessary directories etc. globals may be remade in individual file for compartmentalization

1. rf_cleaning-public cleans the data for use in R and runs logistic/linear reg in stata, some descriptive

2. CF_analysis runs the entire CF analysis pipeline (this can be slow)

3. table_refinement cleans up the R tables into closer to the publication ready versions
-a slight caveat is that the N sizes are manually added to the ATE/ATT CF tables from the log file

NOTE ON PATHS: each of the three scripts sets its own working directory at the
top, so those must be edited to match ${wd} below before running. In
particular CF_analysis.R is NOT driven by the Stata globals: open it and set

    wd <- "C:\\Users\\pardonic\\Documents\\HS employment"

(backslashes doubled, or use forward slashes) to the same folder as ${wd}.
Every other path in CF_analysis.R is relative to that, so nothing else in the
R script needs changing.
*/

**# set up
clear
set more off
cap log close

glo wd "C:/Users/pardonic/Documents/HS employment"  //**************************************set user working directory here

cd "${wd}"

//*******************************************************set R location here, should be in program files
glo rexe "C:/Program Files/R/R-4.4.2/bin/R.exe"

* project sub-directories (globals re-defined inside the individual files too)
glo data   "${wd}/data"               // input .dta
glo temp   "${wd}/temp"               // intermediate files and logs
glo output "${wd}/output"             // tables and figures
glo pub    "${wd}/submission_files"   // publication ready copies

* create anything that is missing; capture so existing folders are left alone
cap mkdir "${temp}"
cap mkdir "${output}"
cap mkdir "${output}/by_cf"
cap mkdir "${output}/f1_cf"
cap mkdir "${pub}"

**# Call files

do "${wd}/rf_cleaning-public.do"

* R CMD BATCH writes its own console transcript to the .Rout file; the script
* additionally sinks to temp/CF_analysis_log.txt and temp/CF_analysis_errors.txt
qui: shell "${rexe}" CMD BATCH --no-save --no-restore "${wd}/CF_analysis.R" "${temp}/CF_analysis.Rout"

do "${wd}/table_refinement.do"
