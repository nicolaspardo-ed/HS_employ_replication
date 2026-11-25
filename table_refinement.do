/*
This file modifies R output/excels and moves or copies
files to a publication ready folder


*/


clear all
set more off
set scheme s1mono

loc wd "C:\Users\pardonic\Documents\HS employment"   // set working directory path here before running
cd "`wd'"

glo pub "`wd'\submission_files"
global output "`wd'/output"

*make the submission directory
capture mkdir "${pub}"

**# copy the exisiting files
copy "${output}/table5_conditional_means_intensity.xlsx" "${pub}/", replace

copy "${output}/table_baseline_binary.csv" "${pub}/regression_binary.csv", replace

copy "${output}/baseline_intensity_results.csv" "${pub}/regression_ordinal.csv", replace

copy "${output}/overlap_combined.tiff" "${pub}/figure1.tiff", replace

copy "${output}/f1_cf/cate_het_Y_bachelors_plus_bytstrel.tiff" "${pub}/figure2.tiff", replace

copy "${output}/f1_cf/cate_het_Y_bachelors_plus_byses2.tiff" "${pub}/figure3.tiff", replace

copy "${output}/f1_cf/cate_het_Y_bachelors_plus_bytxcstd.tiff" "${pub}/figure4.tiff", replace

copy "${output}/f1_cf/cate_het_Y_income_bytstrel.tiff" "${pub}/figure5.tiff", replace

copy "${output}/f1_cf/cate_het_Y_income_byses2.tiff" "${pub}/figure6.tiff", replace

copy "${output}/f1_cf/cate_het_Y_income_byxtracu.tiff" "${pub}/figure7.tiff", replace

copy "${output}/by_cf/results.csv" "${pub}/results_binary10.csv", replace 

copy "${output}/f1_cf/results.csv" "${pub}/results_binary12.csv", replace

copy "${output}/f1_cf/macf_ate.csv" "${pub}/results_macf.csv", replace


**# Re-quote the esttab note rows

*esttab wraps every cell as ="...", which puts the quote in second position, so
*a cell. The data rows hold no commas, so only the note rows are affected.
*Standard CSV quoting keeps each note in one cell with its commas intact.

cap program drop requote_notes
program define requote_notes
	syntax using/

	tempfile tmp
	file open fin  using "`using'", read  text
	file open fout using "`tmp'",   write text replace

	file read fin line
	while r(eof) == 0 {
		* a note row is one ="..." cell on its own line; a data row holds
		* several, so it carries the ",=" separator
		if substr(`"`line'"', 1, 2) == `"=""' & !strpos(`"`line'"', `"",=""') {
			loc txt = substr(`"`line'"', 3, length(`"`line'"') - 3)
			* double any internal quote so the field stays well formed
			loc txt = subinstr(`"`txt'"', `"""', `""""', .)
			file write fout `"""' `"`txt'"' `"""' _n
		}
		else file write fout `"`line'"' _n

		file read fin line
	}

	file close fin
	file close fout
	copy "`tmp'" "`using'", replace
end

requote_notes using "${pub}/regression_binary.csv"
requote_notes using "${pub}/regression_ordinal.csv"

**# Rework R output: balance tables

*small utility program for each file

cap program drop balance_xl
program define balance_xl
	syntax using/, Outfile(string) [Treat(string) Control(string) ///
		Arm(string) Note(string)]

	if "`treat'"   == "" loc treat   "Work"
	if "`control'" == "" loc control "No Work"

	import delimited "`using'", clear varnames(1) case(preserve)

	* multi-arm balance tables stack both comparisons in one file
	if "`arm'" != "" {
		keep if comparison == "`arm'"
	}

	loc keepvars label ///
		mean_treat_base se_treat_base mean_ctrl_base se_ctrl_base smd_base ///
		mean_treat_ipw se_treat_ipw mean_ctrl_ipw se_ctrl_ipw smd_ipw
	keep `keepvars'
	order `keepvars'

	* import delimited reads these as float, rounding to 4 decimals
	foreach v of varlist mean* se* smd* {
		recast double `v'
		qui replace `v' = round(`v', 0.0001)
	}

	* shorten the long dummy labels
	replace label = subinstr(label, "Family composition:", "Family:", 1)

	lab var label           "Covariate"
	lab var mean_treat_base "`treat'"
	lab var se_treat_base   "SE"
	lab var mean_ctrl_base  "`control'"
	lab var se_ctrl_base    "SE"
	lab var smd_base        "SMD"
	lab var mean_treat_ipw  "`treat'"
	lab var se_treat_ipw    "SE"
	lab var mean_ctrl_ipw   "`control'"
	lab var se_ctrl_ipw     "SE"
	lab var smd_ipw         "SMD"

	* row 1 = spanners, row 2 = column headers, row 3 onward = data
	loc first = 3
	loc last  = 2 + _N
	loc noterow = `last' + 1

	export excel using "`outfile'", cell(A2) firstrow(varlabels) replace

	putexcel set "`outfile'", modify
	putexcel (B1:F1) = "Before Weighting", merge hcenter
	putexcel (G1:K1) = "IPW Weighted",     merge hcenter
	putexcel (B`first':K`last'), nformat(number_d2)
	if "`note'" != "" putexcel A`noterow' = "Note. `note'"
	putexcel clear
end


loc balnote "SMD is standardized mean difference. Demographics before weighting incorporate survey weighting. IPW weighting incorporate both survey and inverse-probability weights."

* binary treatments
balance_xl using "${output}/by_cf/balance_table.csv", ///
	outfile("${pub}/10th_work_balance.xlsx") ///
	note("`balnote'")

balance_xl using "${output}/f1_cf/balance_table.csv", ///
	outfile("${pub}/12th_work_balance.xlsx") ///
	note("`balnote'")

* work intensity (multi-arm), one sheet per comparison
balance_xl using "${output}/f1_cf/balance_macf_low.csv", ///
	outfile("${pub}/12th_intensity_balance_low.xlsx") ///
	arm("Arm 1 vs 0") treat("1-20 hrs") control("No Work") ///
	note("`balnote'")

balance_xl using "${output}/f1_cf/balance_macf_high.csv", ///
	outfile("${pub}/12th_intensity_balance_high.xlsx") ///
	arm("Arm 2 vs 0") treat("21+ hrs") control("No Work") ///
	note("`balnote'")


**# Rework R output: BLP tables

*collapses the intercept and slope rows R writes per moderator into one row,
*with the SE in the estimate cell; -blp_prep- reshapes one file, -blp_xl- writes one sheet

cap program drop blp_prep
program define blp_prep
	syntax using/, Saving(string) [SUffix(string)]

	if "`suffix'" == "" loc suffix "low"

	import delimited "`using'", clear varnames(1) case(preserve)

	keep source_variable term Estimate Std_Error p_adj_BH

	* p_adj_BH is "NA" on the intercept rows, so it arrives as a string
	cap confirm string variable p_adj_BH
	if !_rc qui destring p_adj_BH, replace force

	foreach v of varlist Estimate Std_Error {
		recast double `v'
		qui replace `v' = round(`v', 0.001)
	}

	* the slope term is named after the moderator in the binary file and "x"
	* in the multi-arm files, so classify off the intercept instead
	gen byte slope = (term != "(Intercept)")
	drop term

	* estimate (SE) in one cell; %12.3f pads, so trim both sides
	gen cell = trim(string(Estimate, "%12.3f")) + ///
		" (" + trim(string(Std_Error, "%12.3f")) + ")"
	drop Estimate Std_Error

	* R writes the moderators in importance order; keep it
	gen long row = _n
	egen long ord = min(row), by(source_variable)
	drop row

	qui reshape wide cell p_adj_BH, i(source_variable ord) j(slope)

	drop p_adj_BH0
	ren cell0     int_`suffix'
	ren cell1     est_`suffix'
	ren p_adj_BH1 p_`suffix'
	ren ord       ord_`suffix'

	qui save "`saving'", replace
end


cap program drop blp_xl
program define blp_xl
	syntax using/, Outfile(string) SHeet(string) LAbfile(string) ///
		[Highfile(string) SPANLow(string) SPANHigh(string) ///
		Note(string) FIRSTsheet]

	tempfile low high
	blp_prep using "`using'", saving("`low'") suffix(low)

	loc twoarm = ("`highfile'" != "")

	if `twoarm' {
		blp_prep using "`highfile'", saving("`high'") suffix(high)
		use "`low'", clear
		qui merge 1:1 source_variable using "`high'"
		qui count if _merge != 3
		if r(N) > 0 {
			di as err "moderator present in only one arm:"
			list source_variable _merge if _merge != 3
		}
		drop _merge
	}
	else use "`low'", clear

	* covariate labels come from the balance table rather than being typed
	ren source_variable variable
	qui merge m:1 variable using "`labfile'", keep(master match) nogen
	qui replace label = variable if missing(label)

	sort ord_low
	if `twoarm' {
		order label int_low est_low p_low int_high est_high p_high
	}
	else order label int_low est_low p_low
	drop variable ord_*

	lab var label   "Covariate"
	lab var int_low "Intercept"
	lab var est_low "Estimate"
	lab var p_low   "Adjusted p-value"
	if `twoarm' {
		lab var int_high "Intercept"
		lab var est_high "Estimate"
		lab var p_high   "Adjusted p-value"
	}

	* two arms: row 1 = spanners, row 2 = headers, row 3 on = data
	* one arm: no spanner needed, so headers go in row 1
	loc top = cond(`twoarm', 2, 1)
	loc first = `top' + 1
	loc last  = `top' + _N
	loc noterow = `last' + 1

	if "`firstsheet'" != "" {
		export excel using "`outfile'", sheet("`sheet'") ///
			cell(A`top') firstrow(varlabels) replace
	}
	else {
		export excel using "`outfile'", sheet("`sheet'", modify) ///
			cell(A`top') firstrow(varlabels)
	}

	putexcel set "`outfile'", sheet("`sheet'") modify
	if `twoarm' {
		putexcel (B1:D1) = "`spanlow'",  merge hcenter
		putexcel (E1:G1) = "`spanhigh'", merge hcenter
	}
	putexcel (D`first':D`last'), nformat("0.000")
	if `twoarm' putexcel (G`first':G`last'), nformat("0.000")
	if "`note'" != "" putexcel A`noterow' = "note: `note'"
	putexcel clear

	di as txt "wrote `outfile' sheet `sheet' (`=_N' moderators)"
end


* one lookup of covariate labels, shared by all three sheets
tempfile covlab
import delimited "${output}/f1_cf/balance_table.csv", clear varnames(1) case(preserve)
keep variable label
qui duplicates drop variable, force
qui save "`covlab'"

loc blpnote "Standard errors in parentheses. p values are Benjamini and Hochberg adjusted. Intercept is the ATE when the variable is 0. Estimate is the coefficient for a 1 unit increase in the variable."

* binary 12th grade treatment, income (no BA table: RATE was not significant)
blp_xl using "${output}/f1_cf/blp_INC.csv", ///
	outfile("${pub}/blp_tables.xlsx") sheet("12th_income") ///
	labfile("`covlab'") firstsheet note("`blpnote'")

* work intensity, low and high arms side by side
blp_xl using "${output}/f1_cf/blp_macf_BA_low.csv", ///
	highfile("${output}/f1_cf/blp_macf_BA_high.csv") ///
	outfile("${pub}/blp_tables.xlsx") sheet("intensity_BA") ///
	labfile("`covlab'") spanlow("Low Intensity") spanhigh("High Intensity") ///
	note("`blpnote'")

blp_xl using "${output}/f1_cf/blp_macf_INC_low.csv", ///
	highfile("${output}/f1_cf/blp_macf_INC_high.csv") ///
	outfile("${pub}/blp_tables.xlsx") sheet("intensity_INC") ///
	labfile("`covlab'") spanlow("Low Intensity") spanhigh("High Intensity") ///
	note("`blpnote'")



**# rework results

*two-panel ATE/ATT table for the binary treatments: panel 1 = 10th grade,
*panel 2 = 12th grade, rows = outcomes, columns = ATE and ATT

import delimited "${pub}/results_binary10.csv", clear

replace v1 = 10 //grade note to avoid confusion
tempfile r1
save `r1'

import delimited "${pub}/results_binary12.csv", clear
replace v1 = 12

append using `r1'

*calculating p-values for sig stars
gen z_stat_ate = ate/ate_se
gen z_stat_att = att/att_se

gen p_ate = 2 * (1 - normal(abs(z_stat_ate)))

gen p_att = 2 * (1 - normal(abs(z_stat_att)))

gsort v1 -outcome

drop v1 z*

foreach v of varlist ate* att*{
	tostring `v', replace format(%9.3f) force
}

*put in cell together; tostring pads to the %9 width, so trim both pieces
replace ate = strtrim(ate) + " (" + strtrim(ate_se) + ")"
replace att = strtrim(att) + " (" + strtrim(att_se) + ")"

drop *_se

replace ate = ate + "*" if p_ate < 0.05
replace ate = ate + "*" if p_ate < 0.01
replace ate = ate + "*" if p_ate < 0.001

replace att = att + "*" if p_att < 0.05
replace att = att + "*" if p_att < 0.01
replace att = att + "*" if p_att < 0.001


replace outcome = "BA+" if outcome == "BA or Higher"
replace outcome = "Income" if outcome == "Annual Earnings ($)"

drop p_*

*layout below assumes two outcomes per grade, in the order set by -gsort-
assert _N == 4

loc atenote "Standard errors in parentheses. ATE is average treatment effect. ATT is average treatment effect on the treated. * p<0.05, ** p<0.01, *** p<0.001"

*putexcel
putexcel set "${pub}/ate_multipanel.xlsx", replace

putexcel (B1:C1) = "10th Grade", merge hcenter  //header for panel A
putexcel (B5:C5) = "12th Grade", merge hcenter  //header for panel B
putexcel B2 = "ATE"
putexcel C2 = "ATT"

*obs 1-2 go to rows 3-4, obs 3-4 to rows 6-7 (the panel B header takes row 5)
forvalues i = 1/4 {
	loc xr = cond(`i' <= 2, `i' + 2, `i' + 3)
	putexcel A`xr' = (outcome[`i'])
	putexcel B`xr' = (ate[`i'])
	putexcel C`xr' = (att[`i'])
}

putexcel A8 = "Note. `atenote'"
putexcel clear

**# Rework multi-arm results

import delimited "${pub}/results_macf.csv", clear

lab var outcome ""
lab var comparison "Contrast"
lab var ate_base "ATE"
lab var ate_trim "ATE Trimmed"

gen z_stat_ate = ate_base/se_base
gen z_stat_atetrim = ate_trim/se_trim

gen p_atebase = 2 * (1 - normal(abs(z_stat_ate)))

gen p_atetrim = 2 * (1 - normal(abs(z_stat_atetrim)))

drop z*

foreach v of varlist ate* se*{
	tostring `v', replace format(%9.3f) force
	 }

foreach s in base trim{
	replace ate_`s' = ate_`s' + " (" + se_`s' + ")"
	replace ate_`s' = ate_`s' + "*" if p_ate`s' < 0.05
	replace ate_`s' = ate_`s' + "*" if p_ate`s' < 0.01
	replace ate_`s' = ate_`s' + "*" if p_ate`s' < 0.001
}

drop se* p*

replace outcome = "BA+" if outcome == "BA or Higher"
replace outcome = "Income" if outcome == "Annual Earnings ($)"

insobs 1, after(_N)

replace outcome = "Note. Standard errors in parentheses. ATE is average treatment effect. Trimmed sample drops observations with propensities < 0.1. * p<0.05, ** p<0.01, *** p<0.001" if _n == 5

export excel "${pub}/ate_macf_clean.xlsx", firstrow(varlabels) replace 