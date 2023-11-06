**# Notes
*** Characteristics variables
* pseudo_key: anonymized pseudo-id of patients
* date_index: index date, defined as date of SARS-CoV-2 infection diagnosis or symptom onset, whichever occurred earlier
* date_paxlovid: Date of nirmatrelvir/ritonavir initiation
* group: 1 for early nirmatrelvir/ritonavir users (initiated within 1 day from index date); and 0 refers to late nirmatrelvir/ritonavir users (on or beyond 2 days from index date)
* age: Age of patients (years)
* sex: Sex of patients (1 for male; 0 for female)
* charlson_index: Charlson's index of patients
* symptomatic_bl: Symptomatic presentation at baseline
* steroid_covid_bl: Concomitant corticosteroid use at baseline
* immuno_hist: 1 for immunocompromised; 0 for otherwise
* healthcare: Healthcare utilization (inpatient and/or outpatient encounters) in the past year
* covid_hist: Previous SARS-CoV-2 infection
* vaccine_status: COVID-19 vaccination status (1 for "Not fully vaccinated"; 2 for "Fully vaccinated but not boosted"; 3 for "Boosted")
* month_covid_cate: Time of SARS-CoV-2 infection diagnosis (1 for March-June 2022; 2 for July-October 2022; 3 for November 2022-January 2023)
* detection_cate: Type of viral test for case detection (1 for Rapid antigen test; 2 for RT-qPCR)
* dow_covid: Day of the week of SARS-CoV-2 infection diagnosis
* weekend_covid_cate: Day of the week of SARS-CoV-2 infection diagnosis
* ctv_nm: 1 for "With at least one Ct value measurement"; 0 for "Without any Ct value measurements"
* inpatient_paxlovid: 1 for inpatient nirmatrelvir/ritonavir users; 0 for outpatient nirmatrelvir/ritonavir users

*** Study outcomes
* date_death: Date of registered death
* date_admission: Date of hospital admission
* date_VBR: Date of viral burden rebound

*** Datasets
* "covid paxlovid clinical main.dta": list of patients eligible for the mortality or hospitalization analysis
* "covid paxlovid VBR main.dta": list of patients eligible for the viral burden rebound analysis
* "covid paxlovid characteristics.dta": characteristics variables for all analyses

********************************************************************************

**# Table 1. Baseline characteristics of nirmatrelvir/ritonavir users who initiated early (day 0-1) and late (days ≥2) in (a) mortality or hospitalization analysis and (b) viral burden rebound analysis
**# Mortality or hospitalization analysis - IPW
qui forvalues k = 1/18 {
	if `k' == 1 {
		noi di "subgp" _col(15) "gp" _col(30) "mean" _col(45) "sd" _col(60) "min" _col(75) "p25" _col(90) "p50" _col(105) "p90" _col(120) "max"
	}
	use "covid paxlovid clinical main.dta", clear
	keep if subgp_`k' == 1
	merge m:1 pseudo_key using "covid paxlovid characteristics.dta", keep(3) nogen

	misstable sum group age gender charlson_index steroid_covid_bl immuno_hist covid_hist symptomatic_bl healthcare region_num_2 region_num_3 region_num_4 vaccine_status_2 vaccine_status_3 month_covid_cate_2 month_covid_cate_3 detection_cate
	* paxlovid
	logit group age gender charlson_index steroid_covid_bl immuno_hist covid_hist symptomatic_bl healthcare region_num_2 region_num_3 region_num_4 vaccine_status_2 vaccine_status_3 month_covid_cate_2 month_covid_cate_3 detection_cate

	predict pr
	gen pscore = pr*group + (1-pr)*(1-group)

	* generate denominator
	count
	scalar N_all = r(N)
	bysort group : gen N_group = _N
	gen p_i = N_group / N_all
	tab p_i group

	* generate IPTW
	gen _iptw = 1 / pscore
	sum _iptw if group == 1, d
	sum _iptw if group == 0, d

	* generate stabilized IPTW
	gen _iptw_stab = p_i / pscore
	sum _iptw_stab if group == 1, d
	noi di "subgp_`k'" _col(12) "gp=1" _col(20) "N="r(N) _col(30) r(mean) _col(45) r(sd) _col(60) r(min) _col(75) r(p25) _col(90) r(p50) _col(105) r(p90) _col(120) r(max)
	sum _iptw_stab if group == 0, d
	noi di "subgp_`k'" _col(12) "gp=0" _col(20) "N="r(N) _col(30) r(mean) _col(45) r(sd) _col(60) r(min) _col(75) r(p25) _col(90) r(p50) _col(105) r(p90) _col(120) r(max)
	sort pseudo_key
	compress
	save "covid paxlovid clinical IPTW subgp_`k'.dta", replace
}
*
cls
**# Table 1 - mortality or hospitalization analysis - before IPW
qui forvalues k = 1(-1)0 {
	use "covid paxlovid clinical main", clear
	noi di _newline "group=`k'" _newline "var" _col(30) "N" _col(45) "mean" _col(60) "sd"
	foreach var in age charlson_index {
		sum `var' if group == `k'
		noi di "`var'" _col(30) r(N) _col(45) r(mean) _col(60) r(sd)
	}
	noi di "var" _col(30) "N" _col(45) "%"
	qui foreach var in symptomatic_bl steroid_covid_bl immuno_hist healthcare covid_hist {
		count if group == `k'
		scalar N_all = r(N)
		count if `var' == 1 & group == `k'
		scalar N = r(N)
		noi di "`var'=1" _col(30) N _col(45) N/N_all*100
	}
	foreach var in age_gp gender cci_gp vaccine_status month_covid_cate region_num detection_cate dow_covid weekend_covid_cate {
		if inlist("`var'", "age_gp", "cci_gp", "vaccine_status", "month_covid_cate") {
		forvalues j = 1/3 {
			count if group == `k'
			scalar N_all = r(N)
			count if `var' == `j' & group == `k'
			scalar N = r(N)
			noi di "`var'=`j'" _col(30) N _col(45) N/N_all*100
		}
		}
		if inlist("`var'", "region_num") {
		forvalues j = 1/4 {
			gen `var'_`j' = `var' == `j'
			count if group == `k'
			scalar N_all = r(N)
			sum `var'_`j' if group == `k'
			scalar prop = r(mean)
			noi di "`var'=`j'" _col(30) prop*N_all _col(45) prop*100
		}
		}
		if inlist("`var'", "dow_covid") {
		forvalues j = 0/6 {
			gen `var'_`j' = `var' == `j'
			count if group == `k'
			scalar N_all = r(N)
			sum `var'_`j' if group == `k'
			scalar prop = r(mean)
			noi di "`var'=`j'" _col(30) prop*N_all _col(45) prop*100
		}
		}
		if inlist("`var'", "gender", "detection_cate", "weekend_covid_cate") {
		forvalues j = 1(-1)0 {
			count if group == `k'
			scalar N_all = r(N)
			count if `var' == `j' & group == `k'
			scalar N = r(N)
			noi di "`var'=`j'" _col(30) N _col(45) N/N_all*100
		}
		}
	}
}
*
qui forvalues k = 1/1 {
	use "covid paxlovid clinical main", clear
	noi di "var" _col(30) "smd" _col(45) "P"
	foreach var in age charlson_index {
		stddiff `var', by(group) abs cohensd
		scalar SMD = abs(r(stddiff)[1,1])
		regress `var' i.group
		scalar P_value = r(table)[4,2]
		noi di "`var'" _col(30) SMD _col(45) P_value
	}

	foreach var in symptomatic_bl steroid_covid_bl immuno_hist healthcare covid_hist age_gp gender cci_gp vaccine_status month_covid_cate region_num detection_cate dow_covid weekend_covid_cate {
		tab `var' group
		if r(r) > 1 { 
			stddiff i.`var', by(group) abs cohensd
			scalar SMD = abs(r(stddiff)[1,1])
			if ~inlist("`var'", "age_gp", "cci_gp", "vaccine_status", "month_covid_cate", "region_num", "dow_covid") {
				logit `var' i.group
				scalar P_value = r(table)[4,2]
			}
			if inlist("`var'", "age_gp", "cci_gp") {
				ologit `var' i.group
				test 1.group
				scalar P_value = r(p)
			}
			if inlist("`var'", "vaccine_status", "month_covid_cate", "region_num", "modeofcasedetection_num", "dow_covid") {
				capture mlogit `var' i.group
				test 1.group
				scalar P_value = r(p)
			}
			noisily di "`var'" _col(30) SMD _col(45) P_value
		}
		else {
			noisily di "`var'" _col(30) .
		}
	}
}
*
cls
**# Table 1 - mortality or hospitalization analysis - after IPW
qui forvalues k = 1(-1)0 {
	use "covid paxlovid clinical IPTW subgp_1.dta", clear

	noi di _newline "group=`k'" _newline "var" _col(30) "N" _col(45) "mean" _col(60) "sd"
	foreach var in age charlson_index {
		sum `var' if group == `k' [w=_iptw_stab]
		noi di "`var'" _col(30) r(N) _col(45) r(mean) _col(60) r(sd)
	}
	noi di "var" _col(30) "N" _col(45) "%"
	qui foreach var in symptomatic_bl steroid_covid_bl immuno_hist healthcare covid_hist {
		sum `var' if group == `k' [w=_iptw_stab]
		noi di "`var'=1" _col(30) r(N)*r(mean) _col(45) r(mean)*100
	}
	foreach var in age_gp gender cci_gp vaccine_status month_covid_cate region_num detection_cate dow_covid weekend_covid_cate {
		if inlist("`var'", "age_gp", "cci_gp", "vaccine_status", "month_covid_cate") {
		forvalues j = 1/3 {
			gen `var'`j' = `var' == `j'
			sum `var'`j' if group == `k' [w=_iptw_stab]
			noi di "`var'=`j'" _col(30) r(N)*r(mean) _col(45) r(mean)*100
		}
		}
		if inlist("`var'", "region_num") {
		forvalues j = 1/4 {
			gen `var'`j' = `var' == `j'
			sum `var'`j' if group == `k' [w=_iptw_stab]
			noi di "`var'=`j'" _col(30) r(N)*r(mean) _col(45) r(mean)*100
		}
		}
		if inlist("`var'", "dow_covid") {
		forvalues j = 0/6 {
			gen `var'`j' = `var' == `j'
			sum `var'`j' if group == `k' [w=_iptw_stab]
			noi di "`var'=`j'" _col(30) r(N)*r(mean) _col(45) r(mean)*100
		}
		}
		if inlist("`var'", "gender", "detection_cate", "weekend_covid_cate") {
		forvalues j = 1(-1)0 {
			gen `var'`j' = `var' == `j'
			sum `var'`j' if group == `k' [w=_iptw_stab]
			noi di "`var'=`j'" _col(30) r(N)*r(mean) _col(45) r(mean)*100
		}
		}
	}
}
*
qui forvalues k = 1/1 {
	use "covid paxlovid clinical IPTW subgp_1.dta", clear
	noi di "var" _col(30) "smd" _col(45) "P"
	foreach var in age charlson_index {
		sum `var' if group == 1 [w=_iptw_stab]
		local mean_1 = r(mean)
		local sd_1 = r(sd)
		sum `var' if group == 0 [w=_iptw_stab]
		local mean_0 = r(mean)
		local sd_0 = r(sd)
		stddiffi `mean_1' `sd_1' `mean_0' `sd_0'
		scalar SMD = abs(r(std_diff))
		regress `var' i.group
		scalar P_value = r(table)[4,2]
		noi di "`var'" _col(30) SMD _col(45) P_value
	}

	foreach var in symptomatic_bl steroid_covid_bl immuno_hist healthcare covid_hist age_gp gender cci_gp vaccine_status month_covid_cate region_num detection_cate dow_covid weekend_covid_cate {
		tab `var' group
		if r(r) > 1 {
			if r(r) == 2 {
				tab `var' if group == 1 [aw=_iptw_stab], matcell(temp)
				local m11 = int(temp[1,1])
				local m21 = int(temp[2,1])
				tab `var' if group == 0 [aw=_iptw_stab], matcell(temp)
				local m12 = int(temp[1,1])
				local m22 = int(temp[2,1])
				capture stddiffi `m11' `m12' \ `m21' `m22'
			}
			if r(r) == 3 {
				tab `var' if group == 1 [aw=_iptw_stab], matcell(temp)
				local m11 = int(temp[1,1])
				local m21 = int(temp[2,1])
				local m31 = int(temp[3,1])
				tab `var' if group == 0 [aw=_iptw_stab], matcell(temp)
				local m12 = int(temp[1,1])
				local m22 = int(temp[2,1])
				local m32 = int(temp[3,1])
				capture stddiffi `m11' `m12' \ `m21' `m22' \ `m31' `m32'
			}
			if r(r) == 4 {
				tab `var' if group == 1 [aw=_iptw_stab], matcell(temp)
				local m11 = int(temp[1,1])
				local m21 = int(temp[2,1])
				local m31 = int(temp[3,1])
				local m41 = int(temp[4,1])
				tab `var' if group == 0 [aw=_iptw_stab], matcell(temp)
				local m12 = int(temp[1,1])
				local m22 = int(temp[2,1])
				local m32 = int(temp[3,1])
				local m42 = int(temp[4,1])
				capture stddiffi `m11' `m12' \ `m21' `m22' \ `m31' `m32' \ `m41' `m42'
			}
			if r(r) == 7 {
				tab `var' if group == 1 [aw=_iptw_stab], matcell(temp)
				local m11 = int(temp[1,1])
				local m21 = int(temp[2,1])
				local m31 = int(temp[3,1])
				local m41 = int(temp[4,1])
				local m51 = int(temp[5,1])
				local m61 = int(temp[6,1])
				local m71 = int(temp[7,1])
				tab `var' if group == 0 [aw=_iptw_stab], matcell(temp)
				local m12 = int(temp[1,1])
				local m22 = int(temp[2,1])
				local m32 = int(temp[3,1])
				local m42 = int(temp[4,1])
				local m52 = int(temp[5,1])
				local m62 = int(temp[6,1])
				local m72 = int(temp[7,1])
				capture stddiffi `m11' `m12' \ `m21' `m22' \ `m31' `m32' \ `m41' `m42' \ `m51' `m52' \ `m61' `m62' \ `m71' `m72'
			}
			scalar SMD = abs(r(std_diff))
			if ~inlist("`var'", "age_gp", "cci_gp", "vaccine_status", "month_covid_cate", "region_num", "dow_covid") {
				logit `var' i.group [pw=_iptw_stab]
				scalar P_value = r(table)[4,2]
			}
			if inlist("`var'", "age_gp", "cci_gp") {
				ologit `var' i.group [pw=_iptw_stab]
				test 1.group
				scalar P_value = r(p)
			}
			if inlist("`var'", "vaccine_status", "month_covid_cate", "region_num", "modeofcasedetection_num", "dow_covid") {
				capture mlogit `var' i.group [pw=_iptw_stab]
				test 1.group
				scalar P_value = r(p)
			}
			noisily di "`var'" _col(30) SMD _col(45) P_value
		}
		else {
			noisily di "`var'" _col(30) .
		}
	}
}
*

**# Viral burden rebound analysis - IPW
qui forvalues k = 1/18 {
	use "covid paxlovid VBR main.dta", clear
	merge m:1 pseudo_key using "covid paxlovid characteristics.dta", keep(3) nogen
	keep if subgp_`k' == 1

	logit group age gender charlson_index steroid_covid_bl immuno_hist covid_hist symptomatic_bl healthcare region_num_2 region_num_3 region_num_4 vaccine_status_2 vaccine_status_3 month_covid_cate_2 month_covid_cate_3 detection_cate inpatient_paxlovid censor_ctv_day_5

	predict prob_tx if e(sample)
	gen pscore_denom = prob_tx*group + (1-prob_tx)*(1-group)

	* generate numerator
	count
	scalar N_all = r(N)
	bysort group : gen N_group = _N
	gen pscore_num = N_group / N_all
	tab pscore_num group

	* generate stabilized IPTW
	gen _iptw_stab = pscore_num / pscore_denom
	sum _iptw_stab if group == 1, d
	noi di "subgp_`k'" _col(12) "gp=1" _col(20) "N="r(N) _col(30) r(mean) _col(45) r(sd) _col(60) r(min) _col(75) r(p25) _col(90) r(p50) _col(105) r(p90) _col(120) r(max)
	sum _iptw_stab if group == 0, d
	noi di "subgp_`k'" _col(12) "gp=0" _col(20) "N="r(N) _col(30) r(mean) _col(45) r(sd) _col(60) r(min) _col(75) r(p25) _col(90) r(p50) _col(105) r(p90) _col(120) r(max)
	sort pseudo_key
	compress
	save "covid paxlovid VBR IPTW subgp_`k'.dta", replace
}
*
cls
**# Table 1 - viral burden rebound analysis - before IPW
qui forvalues k = 1(-1)0 {
use "covid paxlovid VBR main.dta", clear
noi di _newline "group=`k'" _newline "var" _col(30) "N" _col(45) "mean" _col(60) "sd"
foreach var in age charlson_index {
	sum `var' if group == `k'
	noi di "`var'" _col(30) r(N) _col(45) r(mean) _col(60) r(sd)
}
noi di "var" _col(30) "N" _col(45) "%"
qui foreach var in symptomatic_bl inpatient_paxlovid steroid_covid_bl immuno_hist healthcare covid_hist censor_ctv_day_5_0 {
	count if group == `k'
	scalar N_all = r(N)
	count if `var' == 1 & group == `k'
	scalar N = r(N)
	noi di "`var'=1" _col(30) N _col(45) N/N_all*100
}
foreach var in age_gp gender cci_gp vaccine_status month_covid_cate region_num detection_cate dow_covid weekend_covid_cate {
	if inlist("`var'", "age_gp", "cci_gp", "vaccine_status", "month_covid_cate") {
	forvalues j = 1/3 {
		count if group == `k'
		scalar N_all = r(N)
		count if `var' == `j' & group == `k'
		scalar N = r(N)
		noi di "`var'=`j'" _col(30) N _col(45) N/N_all*100
	}
	}
	if inlist("`var'", "region_num") {
	forvalues j = 1/4 {
		gen `var'_`j' = `var' == `j'
		count if group == `k'
		scalar N_all = r(N)
		sum `var'_`j' if group == `k'
		scalar prop = r(mean)
		noi di "`var'=`j'" _col(30) prop*N_all _col(45) prop*100
	}
	}
	if inlist("`var'", "dow_covid") {
	forvalues j = 0/6 {
		gen `var'_`j' = `var' == `j'
		count if group == `k'
		scalar N_all = r(N)
		sum `var'_`j' if group == `k'
		scalar prop = r(mean)
		noi di "`var'=`j'" _col(30) prop*N_all _col(45) prop*100
	}
	}
	if inlist("`var'", "gender", "detection_cate", "weekend_covid_cate") {
	forvalues j = 1(-1)0 {
		count if group == `k'
		scalar N_all = r(N)
		count if `var' == `j' & group == `k'
		scalar N = r(N)
		noi di "`var'=`j'" _col(30) N _col(45) N/N_all*100
	}
	}
}
}
*
qui forvalues k = 1/1 {
use "covid paxlovid VBR main.dta", clear
noi di "var" _col(30) "smd" _col(45) "P"
foreach var in age charlson_index {
	stddiff `var', by(group) abs cohensd
	scalar SMD = abs(r(stddiff)[1,1])
	regress `var' i.group
	scalar P_value = r(table)[4,2]
	noi di "`var'" _col(30) SMD _col(45) P_value
}

foreach var in symptomatic_bl inpatient_paxlovid steroid_covid_bl immuno_hist healthcare covid_hist censor_ctv_day_5 age_gp gender cci_gp vaccine_status month_covid_cate region_num detection_cate dow_covid weekend_covid_cate {
	tab `var' group
	if r(r) > 1 { 
		stddiff i.`var', by(group) abs cohensd
		scalar SMD = abs(r(stddiff)[1,1])
		if ~inlist("`var'", "age_gp", "cci_gp", "vaccine_status", "month_covid_cate", "region_num", "inpatient_paxlovid", "dow_covid") {
			logit `var' i.group
			scalar P_value = r(table)[4,2]
		}
		if inlist("`var'", "age_gp", "cci_gp") {
			ologit `var' i.group
			test 1.group
			scalar P_value = r(p)
		}
		if inlist("`var'", "vaccine_status", "month_covid_cate", "region_num", "detection_cate", "dow_covid") {
			capture mlogit `var' i.group
			test 1.group
			scalar P_value = r(p)
		}
		noisily di "`var'" _col(30) SMD _col(45) P_value
	}
	else {
		noisily di "`var'" _col(30) .
	}
}
}
*
cls
**# Table 1 - viral burden rebound analysis - after IPW
qui forvalues k = 1(-1)0 {
use "covid paxlovid VBR IPTW subgp_1.dta", clear
noi di _newline "group=`k'" _newline "var" _col(30) "N" _col(45) "mean" _col(60) "sd"
foreach var in age charlson_index {
	sum `var' if group == `k' [w=_iptw_stab]
	noi di "`var'" _col(30) r(N) _col(45) r(mean) _col(60) r(sd)
}
noi di "var" _col(30) "N" _col(45) "%"
qui foreach var in symptomatic_bl inpatient_paxlovid steroid_covid_bl immuno_hist healthcare covid_hist censor_ctv_day_5_0 {
	sum `var' if group == `k' [w=_iptw_stab]
	noi di "`var'=1" _col(30) r(N)*r(mean) _col(45) r(mean)*100
}
foreach var in age_gp gender cci_gp vaccine_status month_covid_cate region_num detection_cate dow_covid weekend_covid_cate {
	if inlist("`var'", "age_gp", "cci_gp", "vaccine_status", "month_covid_cate") {
	forvalues j = 1/3 {
		gen `var'`j' = `var' == `j'
		sum `var'`j' if group == `k' [w=_iptw_stab]
		noi di "`var'=`j'" _col(30) r(N)*r(mean) _col(45) r(mean)*100
	}
	}
	if inlist("`var'", "region_num") {
	forvalues j = 1/4 {
		gen `var'`j' = `var' == `j'
		sum `var'`j' if group == `k' [w=_iptw_stab]
		noi di "`var'=`j'" _col(30) r(N)*r(mean) _col(45) r(mean)*100
	}
	}
	if inlist("`var'", "dow_covid") {
	forvalues j = 0/6 {
		gen `var'`j' = `var' == `j'
		sum `var'`j' if group == `k' [w=_iptw_stab]
		noi di "`var'=`j'" _col(30) r(N)*r(mean) _col(45) r(mean)*100
	}
	}
	if inlist("`var'", "gender", "detection_cate", "weekend_covid_cate") {
	forvalues j = 1(-1)0 {
		gen `var'`j' = `var' == `j'
		sum `var'`j' if group == `k' [w=_iptw_stab]
		noi di "`var'=`j'" _col(30) r(N)*r(mean) _col(45) r(mean)*100
	}
	}
}
}
*
qui forvalues k = 1/1 {
use "covid paxlovid VBR IPTW subgp_1.dta", clear
noi di "var" _col(30) "smd" _col(45) "P"
foreach var in age charlson_index {
	sum `var' if group == 1 [w=_iptw_stab]
	local mean_1 = r(mean)
	local sd_1 = r(sd)
	sum `var' if group == 0 [w=_iptw_stab]
	local mean_0 = r(mean)
	local sd_0 = r(sd)
	stddiffi `mean_1' `sd_1' `mean_0' `sd_0'
	scalar SMD = abs(r(std_diff))
	regress `var' i.group
	scalar P_value = r(table)[4,2]
	noi di "`var'" _col(30) SMD _col(45) P_value
}

foreach var in symptomatic_bl inpatient_paxlovid steroid_covid_bl immuno_hist healthcare covid_hist censor_ctv_day_5 age_gp gender cci_gp vaccine_status month_covid_cate region_num detection_cate dow_covid weekend_covid_cate {
	tab `var' group
	if r(r) > 1 {
		if r(r) == 2 {
			tab `var' if group == 1 [aw=_iptw_stab], matcell(temp)
			local m11 = int(temp[1,1])
			local m21 = int(temp[2,1])
			tab `var' if group == 0 [aw=_iptw_stab], matcell(temp)
			local m12 = int(temp[1,1])
			local m22 = int(temp[2,1])
			capture stddiffi `m11' `m12' \ `m21' `m22'
		}
		if r(r) == 3 {
			tab `var' if group == 1 [aw=_iptw_stab], matcell(temp)
			local m11 = int(temp[1,1])
			local m21 = int(temp[2,1])
			local m31 = int(temp[3,1])
			tab `var' if group == 0 [aw=_iptw_stab], matcell(temp)
			local m12 = int(temp[1,1])
			local m22 = int(temp[2,1])
			local m32 = int(temp[3,1])
			capture stddiffi `m11' `m12' \ `m21' `m22' \ `m31' `m32'
		}
		if r(r) == 4 {
			tab `var' if group == 1 [aw=_iptw_stab], matcell(temp)
			local m11 = int(temp[1,1])
			local m21 = int(temp[2,1])
			local m31 = int(temp[3,1])
			local m41 = int(temp[4,1])
			tab `var' if group == 0 [aw=_iptw_stab], matcell(temp)
			local m12 = int(temp[1,1])
			local m22 = int(temp[2,1])
			local m32 = int(temp[3,1])
			local m42 = int(temp[4,1])
			capture stddiffi `m11' `m12' \ `m21' `m22' \ `m31' `m32' \ `m41' `m42'
		}
		if r(r) == 7 {
			tab `var' if group == 1 [aw=_iptw_stab], matcell(temp)
			local m11 = int(temp[1,1])
			local m21 = int(temp[2,1])
			local m31 = int(temp[3,1])
			local m41 = int(temp[4,1])
			local m51 = int(temp[5,1])
			local m61 = int(temp[6,1])
			local m71 = int(temp[7,1])
			tab `var' if group == 0 [aw=_iptw_stab], matcell(temp)
			local m12 = int(temp[1,1])
			local m22 = int(temp[2,1])
			local m32 = int(temp[3,1])
			local m42 = int(temp[4,1])
			local m52 = int(temp[5,1])
			local m62 = int(temp[6,1])
			local m72 = int(temp[7,1])
			capture stddiffi `m11' `m12' \ `m21' `m22' \ `m31' `m32' \ `m41' `m42' \ `m51' `m52' \ `m61' `m62' \ `m71' `m72'
		}
		scalar SMD = abs(r(std_diff))
		
		if ~inlist("`var'", "age_gp", "cci_gp", "vaccine_status", "month_covid_cate", "detection_cate", "region_num", "dow_covid") {
			logit `var' i.group [pw=_iptw_stab]
			scalar P_value = r(table)[4,2]
		}
		if inlist("`var'", "age_gp", "cci_gp") {
			ologit `var' i.group [pw=_iptw_stab]
			test 1.group
			scalar P_value = r(p)
		}
		if inlist("`var'", "vaccine_status", "month_covid_cate", "region_num", "detection_cate", "weekend_covid_cate") {
			capture mlogit `var' i.group [pw=_iptw_stab]
			test 1.group
			scalar P_value = r(p)
		}
		noisily di "`var'" _col(30) SMD _col(45) P_value
	}
	else {
		noisily di "`var'" _col(30) .
	}
}
}
*

************************************************************************************************
**# Table 2. Association between timing of nirmatrelvir/ritonavir initiation and 28-day all-cause mortality or hospitalization.
**# Mortality or hospitalization analysis - Main analysis
* Target Trial Emulation
qui forvalues k = 1/18 {
forvalues j = 0/100 {
	
	*** Prepare dataset for analysis
	capture mkdir "clinical"
	capture mkdir "clinical/subgp_`k'"

	* Setup for trial emulation
	capture mkdir "clinical/subgp_`k'/prepare"
	noi di "subgp_`k'" _col(15) "bs`j'" _col(30) "setup"
	if !fileexists("clinical/subgp_`k'/prepare/prepare subgp_`k' bs`j'.dta") {
	use "covid paxlovid clinical main.dta", clear
	
	* generate subgroups
	gen subgp_1 = 1

	gen subgp_2 = gender == 1
	gen subgp_3 = gender == 0

	gen subgp_4 = month_covid_cate == 1
	gen subgp_5 = month_covid_cate == 2
	gen subgp_6 = month_covid_cate == 3

	gen subgp_7 = vaccine_status > 1
	gen subgp_8 = vaccine_status == 1

	gen subgp_9 = inrange(charlson_index, 0, 6)
	gen subgp_10 = charlson_index > 7

	gen subgp_11 = date_steroid_covid <= date_index
	gen subgp_12 = !(date_steroid_covid <= date_index)

	gen subgp_13 = immuno_hist == 1
	gen subgp_14 = immuno_hist == 0

	gen subgp_15 = ctv_nm == 1
	gen subgp_16 = ctv_nm == 0

	gen subgp_17 = date_index == date_onset
	gen subgp_18 = date_index == date_covid
	
	keep if subgp_`k' == 1
	sort pseudo_key
	if `j' > 0 {
		set seed `j'
		bsample
	}

	gen date_event = min(date_death, date_admission)
	gen time_to_event = date_event - date_index
	gen treatment = inrange(time_to_treatment, 0, 1)
	
	gen censor = 1 if time_to_treatment > 5
	replace censor = 1 if date_paxlovid > date_admission
	
	gen date_censor = date_index + 5 if time_to_treatment > 5
	replace date_censor = min(date_censor, date_admission) if date_paxlovid > date_admission
	
	gen date_start_fu = date_index
	gen date_last_fu = min(date_death, mdy(02,12,2023), date_index + 28, date_molnupiravir, date_censor)
	
	gen event = inrange(date_event, date_start_fu, date_last_fu)
	gen fup_obs = min(date_event-date_index, 28) if event == 1
	replace fup_obs = min(date_last_fu-date_index, 28) if event == 0
	
	format date* %td
	* keep necessary variables
	keep pseudo_key group fup_obs event time_to_treatment time_to_event treatment censor
	gen bs = `j'
	compress
	save "clinical/subgp_`k'/prepare/prepare subgp_`k' bs`j'.dta", replace
	}
	
	*** Cloning & censoring
	noi di "subgp_`k'" _col(15) "bs`j'" _col(30) "cloned" 
	capture mkdir "clinical/subgp_`k'/cloned"
	if !fileexists("clinical/subgp_`k'/cloned/cloned subgp_`k' bs`j'.dta") {
	* Prepare dataset for analysis
	use "clinical/subgp_`k'/prepare/prepare subgp_`k' bs`j'.dta", replace
	stset fup_obs, failure(event)

	* Arm A: no treatment within grace period (control: non-exposed group)
	gen outcomeA = _d // _d = event
	gen fupA = _t // _t = follow up time

	/// if the patient received treatment within day 0-1:
	/// 1. no event outcome, since the patient survived till censoring (treatment)
	replace outcomeA = 0 if treatment == 1 & time_to_treatment <= 1
	/// 2. follow up is censored at treatment
	replace fupA = time_to_treatment if treatment == 1 & time_to_treatment <= 1

	* Arm B: treatment within day 0-1 (treated: exposed group)
	gen outcomeB = _d 
	gen fupB = _t 

	/// if the patient survived the grace period and did not receive treatment within day 0-1:
	/// 1. no event outcome if the patient survived the first day
	replace outcomeB = 0 if (treatment==0 & _t>1) | (treatment==1 & time_to_treatment >1 & time_to_treatment !=.)
	/// 2. follow up is censored at day 1
	replace fupB = 1 if (treatment==0 & _t>1) | (treatment==1 & time_to_treatment >1 & time_to_treatment != .)

	** append clones
	preserve
		drop outcomeB fupB
		rename outcomeA outcome
		rename fupA fup
		gen arm = "NoTreatment"
		tempfile a
		save "`a'", replace
	restore
		drop outcomeA fupA
		rename outcomeB outcome
		rename fupB fup
		gen arm = "Treatment"	
		cap append using "`a'"

	// Weight models

	sort _all
	gen NewID = _n

	** exclude those initiated at day 0 in Late arm
	drop if time_to_treatment == 0 & arm == "NoTreatment"

	***** Weight model: define survival time and event indicator	
	*** Early arm
	** Case 1: they do not deviate at time of treatment, but are not at risk of deviating any more
	gen wm_fup = time_to_treatment if arm == "Treatment" & time_to_treatment<=1 & time_to_treatment!=. & treatment == 1 
	gen wm_outcome = 0 if arm == "Treatment" & time_to_treatment<=1 & time_to_treatment!=. & treatment == 1 

	** Case 2: they deviate at day 1
	replace wm_fup = 1 if arm == "Treatment" & ((treatment == 0 & fup >= 1) | (time_to_treatment>1 & treatment == 1))
	replace wm_outcome = 1 if arm == "Treatment" & ((treatment == 0 & fup >= 1) | (time_to_treatment>1 & treatment == 1)) 

	** Case 3: they do not deviate, but we need to keep their survival as observed and censor them as we do not know what happens afterwards
	replace wm_fup = fup if arm == "Treatment" & treatment == 0 & fup < 1
	replace wm_outcome = 0 if arm == "Treatment" & treatment == 0 & fup < 1
	
	** add 1 days to 0-survivors
	replace wm_fup = 1 if arm == "Treatment" & wm_fup==0 

	*** Late arm
	** Case 1: they do not deviate at time of treatment, but are not at risk of deviating any more
	replace wm_fup = time_to_treatment if arm == "NoTreatment" & time_to_treatment<=1 & treatment == 1 
	replace wm_outcome = 1 if arm == "NoTreatment" & time_to_treatment<=1 & treatment == 1 

	** Case 2: they deviate at day 1
	replace wm_fup = 1 if arm == "NoTreatment" & ((treatment == 0 & fup >= 1) | (time_to_treatment>1 & treatment == 1)) 
	replace wm_outcome = 0 if arm == "NoTreatment" & ((treatment == 0 & fup >= 1) | (time_to_treatment>1 & treatment == 1)) 

	** Case 3: they do not deviate, but we need to keep their survival as observed and censor them as we do not know what happens afterwards
	replace wm_fup = fup if arm == "NoTreatment" & treatment == 0 & fup < 1
	replace wm_outcome = 0 if arm == "NoTreatment" & treatment == 0 & fup < 1
	
	** add 1 days to 0-survivors
	replace wm_fup = 1 if arm == "NoTreatment" & wm_fup==0
	
	*** Censoring criteria - use tab to check the number of patients in each case and their outcome/censoring assignments
	* Case 1: patients admitted and initiated treatment on day 1: outcome in early, censored in late
	*tab outcome arm if time_to_treatment == 1 & time_to_event == 1
	*tab wm_outcome arm if time_to_treatment == 1 & time_to_event == 1

	* Case 2: patients admitted on day 1 with initiation after admission: outcome in late, censored in early
	*tab outcome arm if time_to_event == 1 & time_to_treatment > time_to_event 
	*tab wm_outcome arm if time_to_event == 1 & time_to_treatment > time_to_event
	replace outcome = 0 if arm == "Treatment" & time_to_event == 1 & time_to_treatment > time_to_event
	
	* Case 3: patients admitted on days 2-5 with initiation on the same day: outcome in late, censored in early
	*tab outcome arm if inrange(time_to_treatment,2,5) & time_to_treatment == time_to_event
	*tab fup arm if inrange(time_to_treatment,2,5) & time_to_treatment == time_to_event
	
	* Case 4: patients intiated treatments after 5 days: censor them at day 5
	* check all time_to_treatment > 5 censored on day 5 (or earlier)
	*tab time_to_treatment time_to_event if inrange(time_to_treatment, 0, 5) & inrange(time_to_event, 0, 5)
	*tab time_to_treatment time_to_event if inrange(time_to_event, 0, 5)
	*tab fup arm if time_to_treatment > 5 & time_to_event > 5
	*tab time_to_treatment if time_to_event <= 5
	*tab outcome arm if time_to_treatment > 5 & time_to_event > 5
	*tab wm_outcome arm if time_to_treatment > 5 & time_to_event > 5
	replace outcome = 0 if time_to_treatment > 5 & time_to_event > 5
	replace wm_outcome = 1 if time_to_treatment > 5 & time_to_event > 5
	
	gen arm_value = 1 if arm == "Treatment"
	replace arm_value = 0 if arm == "NoTreatment"
	drop arm
	order pseudo_key Early arm_value
	compress
	save "clinical/subgp_`k'/cloned/cloned subgp_`k' bs`j'.dta", replace
	}
	
	*** Split times
	capture mkdir "clinical/subgp_`k'/split"
	noi di "subgp_`k'" _col(15) "bs`j'" _col(30) "split"
	if !fileexists("clinical/subgp_`k'/split/split subgp_`k' bs`j'.dta") {
	use "clinical/subgp_`k'/cloned/cloned subgp_`k' bs`j'.dta", clear
	
	** times of event
	stset fup, fail(outcome) id(NewID)
	stsplit, at(failures)

	gen tstart = _t0
	gen tstop = _t

	** times of censoring
	gen TrialEmul_cens = 1-outcome
	stset fup, fail(TrialEmul_cens) id(NewID)
	stsplit, at(failures)

	replace tstart = _t0 if tstart<_t0 & _t0 != . & _t != .
	replace tstop = _t if tstop>_t & _t0 != . & _t != .

	order pseudo_key arm tstart tstop
	sort _all
	compress
	save "clinical/subgp_`k'/split/split subgp_`k' bs`j'.dta", replace
	}

	*** IPCW - Early arm
	capture mkdir "clinical/subgp_`k'/Treatment"
	noi di "subgp_`k'" _col(15) "bs`j'" _col(30) "Treatment" 
	if !fileexists("clinical/subgp_`k'/Treatment/Treatment subgp_`k' bs`j'.dta") {
	use "clinical/subgp_`k'/split/split subgp_`k' bs`j'.dta", clear
	keep if arm_value == 1
	merge m:1 pseudo_key using "covid paxlovid characteristics.dta", keep(3) nogen
	replace outcome = . if TrialEmul_cens == .

	* adapt records to the long format
	sort NewID tstart
	bysort NewID (tstart): replace wm_outcome = 0 if _n!=_N
	
	* Weight model:
	stset tstop, origin(time tstart) failure(wm_outcome) id(NewID)
	stcox age sex charlson_index steroid_covid_bl immuno_hist covid_hist symptomatic_bl healthcare region_num_2 region_num_3 region_num_4 vaccine_status_2 vaccine_status_3 month_covid_cate_2 month_covid_cate_3 detection_cate, efron
	if _rc == 0 {
		predict ch2_2, basech
		gen s0_2=exp(-ch2_2)
		predict xbCox2_2, xb
		gen weight = 1/s0_2^(exp(xbCox2_2)) 
		drop _st-xbCox2_2 
		gen invalid = 0
	}
	else {
		gen invalid = 1
		gen weight = .
		noi di "invalid"
	}
	keep pseudo_key arm_value tstart tstop treatment bs outcome wm_outcome NewID weight invalid
	compress
	save "clinical/subgp_`k'/Treatment/Treatment subgp_`k' bs`j'.dta", replace
	}
	
	*** IPCW - Late arm
	capture mkdir "clinical/subgp_`k'/NoTreatment"
	noi di "subgp_`k'" _col(15) "bs`j'" _col(30) "NoTreatment" 
	if !fileexists("clinical/subgp_`k'/NoTreatment/Notreatment subgp_`k' bs`j'.dta") {
	use "clinical/subgp_`k'/split/split subgp_`k' bs`j'.dta", clear
	keep if arm_value == 0
	merge m:1 pseudo_key using "covid paxlovid characteristics.dta", keep(3) nogen
	replace outcome = . if TrialEmul_cens == .

	* adapt records to the long format
	sort NewID tstart
	bysort NewID: replace wm_outcome = 0 if _n!=_N

	* Weight model:
	stset tstop, origin(time tstart) failure(wm_outcome) id(NewID)
	stcox age sex charlson_index steroid_covid_bl immuno_hist covid_hist symptomatic_bl healthcare region_num_2 region_num_3 region_num_4 vaccine_status_2 vaccine_status_3 month_covid_cate_2 month_covid_cate_3 detection_cate, efron
	if _rc == 0 {
		predict ch2_2, basech
		gen s0_2=exp(-ch2_2)
		predict xbCox2_2, xb
		gen weight = 1/s0_2^(exp(xbCox2_2)) 
		drop _st-xbCox2_2 
		gen invalid = 0
	}
	else {
		gen invalid = 1
		gen weight = .
		noi di "invalid"
	}
	keep pseudo_key arm_value tstart tstop treatment bs outcome wm_outcome NewID weight invalid
	compress
	save "clinical/subgp_`k'/NoTreatment/Notreatment subgp_`k' bs`j'.dta", replace
	}
	
	*** Combine & Generate weights
	capture mkdir "clinical/subgp_`k'/combined"
	noi di "subgp_`k'" _col(15) "bs`j'" _col(30) "Combine & weights"
	if !fileexists("clinical/subgp_`k'/weight/weight subgp_`k' bs`j'.dta") {
	use "clinical/subgp_`k'/Treatment/Treatment subgp_`k' bs`j'.dta", clear
	append using "clinical/subgp_`k'/NoTreatment/Notreatment subgp_`k' bs`j'.dta"

	// create a new ID variable for each clone in each arm
	tostring NewID, gen(NewID_str)
	replace NewID_str = "00000" + NewID_str if length(NewID_str)== 1
	replace NewID_str = "0000" + NewID_str if length(NewID_str)== 2
	replace NewID_str = "000" + NewID_str if length(NewID_str)== 3
	replace NewID_str = "00" + NewID_str if length(NewID_str)== 4
	replace NewID_str = "0" + NewID_str if length(NewID_str)== 5
	gen Anal_ID = "1" + NewID_str if arm_value == 1
	replace Anal_ID = "2" + NewID_str if arm_value == 0

	replace weight = 1 if wm_outcome == 1 & tstop == 1
	sum weight
	local max_weight = r(max)

	keep pseudo_key arm_value tstart tstop treatment bs outcome wm_outcome Anal_ID weight invalid
	destring Anal_ID, replace

	* generate denominator
	rename weight _ipcw
	gen prob_denom = 1/_ipcw
	* generate numerator
	bysort arm_value tstart : gen N_group = _N
	bysort arm_value tstart : egen N_censored = sum(wm_outcome)
	gen prob_uncensored = (N_group-N_censored)/N_group
	gen prob_num = prob_uncensored
	sort arm_value pseudo_key tstart
	by arm_value pseudo_key : replace prob_num=prob_num*prob_num[_n-1] if _n!=1
	* generate stabilized IPCW
	gen _ipcw_stab = prob_num / prob_denom
	* generate combined weight
	merge m:1 pseudo_key using "covid paxlovid clinical IPTW subgp_`k'.dta", keep(3) keepusing(_iptw _iptw_stab)
	gen weight = _iptw * _ipcw
	gen weight_stab = _iptw_stab * _ipcw_stab
	compress
	save "clinical/subgp_`k'/weight/weight subgp_`k' bs`j'.dta", replace
	}
	
	*** Generate KM estimate using weight_stab
	if !fileexists("clinical/subgp_`k'/KM weight_stab/KM weight_stab subgp_`k' bs`j'.dta") {
	capture mkdir "clinical/subgp_`k'/KM weight_stab"
	noi di "subgp_`k'" _col(15) "bs`j'" _col(30) "KM weight_stab"
	use "clinical/subgp_`k'/weight/weight subgp_`k' bs`j'.dta", clear
	count if invalid == 1
	if r(N) == 0 {
		stset tstop [pweight = weight_stab], enter(time tstart) failure(outcome)
		sts generate KM_s_w = s if arm_value == 1
		sts generate KM_ns_w = s if arm_value == 0
	}
	else {
		gen KM_s_w = .
		gen KM_ns_w = .
	}
	collapse (firstnm) KM_s_w KM_ns_w, by(tstop bs invalid)
	save "clinical/subgp_`k'/KM weight_stab/KM weight_stab subgp_`k' bs`j'.dta", replace
	}
}
}
*
* Finalize bootstrap datasets
qui forvalues k = 1/18 {
clear
forvalues j = 0/100 {
	capture append using "clinical/subgp_`k'/KM weight_stab/KM weight_stab subgp_`k' bs`j'.dta"
	*erase "unmatched/KM/KM subgp_`k' bs`j'.dta"
}
	gen hazard_s_w = 1 - KM_s_w
	gen hazard_ns_w = 1 - KM_ns_w
	gen odds_s_w = hazard_s_w/(1-hazard_s_w)
	gen odds_ns_w = hazard_ns_w/(1-hazard_ns_w)
	gen RR_w = hazard_s_w/hazard_ns_w
	gen diff_w = hazard_s_w - hazard_ns_w
	gen OR_w = odds_s_w / odds_ns_w
	compress
	save "clinical/subgp_`k'/KM weight_stab subgp_`k' bs_all.dta", replace
}
*
cls
* KM estimate
qui forvalues k = 1/18 {
	use "clinical/subgp_`k'/KM weight_stab subgp_`k' bs_all.dta", clear
	
	keep if KM_s_w < . & KM_ns_w < .
	drop if invalid == 1
	bysort bs (tstop) : keep if _n == _N
	keep if _n <= 1 + 100
	
	sum hazard_s_w if bs == 0, d
	scalar hazard_s_mean = r(mean)
	centile hazard_s_w if bs > 0, centile(2.5 97.5)
	scalar hazard_s_cil = r(c_1)
	scalar hazard_s_ciu = r(c_2)
	
	sum hazard_ns_w if bs == 0, d
	scalar hazard_ns_mean = r(mean)
	centile hazard_ns_w if bs > 0, centile(2.5 97.5)
	scalar hazard_ns_cil = r(c_1)
	scalar hazard_ns_ciu = r(c_2)
	
	noi di "subgp_`k'" _col(15) hazard_s_mean _col(30) hazard_s_cil _col(45) hazard_s_ciu _col(60) ///
	hazard_ns_mean _col(75) hazard_s_cil _col(90) hazard_s_ciu
}
*
* Absolute risk reduction
qui forvalues k = 1/18 {
	use "clinical/subgp_`k'/KM weight_stab subgp_`k' bs_all.dta", clear
	
	keep if KM_s_w < . & KM_ns_w < .
	drop if invalid == 1
	bysort bs (tstop) : keep if _n == _N
	keep if _n <= 1 + 100
	
	sum diff_w if bs == 0, d
	scalar bs_mean = r(mean)
	scalar bs_p50 = r(p50)
	
	centile diff_w if bs > 0, centile(2.5 97.5)
	scalar bs_cil = r(c_1)
	scalar bs_ciu = r(c_2)
	noi di "subgp_`k'" _col(15) bs_mean _col(30) bs_p50 _col(45) bs_cil _col(60) bs_ciu
}
*
* Relative risk
qui forvalues k = 1/18 {
	use "clinical/subgp_`k'/KM weight_stab subgp_`k' bs_all.dta", clear
	
	keep if KM_s_w < . & KM_ns_w < .
	drop if invalid == 1
	bysort bs (tstop) : keep if _n == _N
	keep if _n <= 1 + 100
	
	sum RR_w if bs == 0, d
	scalar bs_mean = r(mean)
	scalar bs_p50 = r(p50)
	centile RR_w if bs > 0, centile(2.5 97.5)
	scalar bs_cil = r(c_1)
	scalar bs_ciu = r(c_2)
	noi di "subgp_`k'" _col(15) bs_mean _col(30) bs_p50 _col(45) bs_cil _col(60) bs_ciu
}
*
* N
qui forvalues k = 1/18 {
	use "clinical/subgp_`k'/cloned/cloned subgp_`k' bs0.dta", clear
	merge m:1 pseudo_key using "covid paxlovid clinical main.dta", keepusing(group) keep(1 3) nogen
	keep if (arm_value == 1 & group == 1) | (arm_value == 0 & group == 0)
	count if group == 1
	scalar N_1 = r(N)
	count if group == 0
	scalar N_0 = r(N)
	count if group == 1 & outcome == 1
	scalar n_e_1 = r(N)
	count if group == 0 & outcome == 1
	scalar n_e_0 = r(N)
	noi di "subgp_`k'" _col(15) N_1 _col(30) N_0 _col(45) n_e_1 _col(60) n_e_0
}
*
qui forvalues k = 1/18 {
	if `k' == 1 {
		noi di "subgp" _col(15) "N_1" _col(30) "N_0" _col(45) "n_e_1" _col(60) "n_e_0"
	}
	use "clinical/subgp_`k'/cloned/cloned subgp_`k' bs0.dta", clear
	merge m:1 pseudo_key using "covid paxlovid clinical IPTW subgp_`k'.dta", keep(3)
	count if arm_value == 1
	scalar N_1 = r(N)
	count if arm_value == 1 & outcome == 1
	scalar n_e_1 = r(N)
	count if arm_value == 0
	scalar N_0 = r(N)
	count if arm_value == 0 & outcome == 1
	scalar n_e_0 = r(N)
	
	noi di "subgp_`k'" _col(15) N_1 _col(30) N_0 _col(45) n_e_1 _col(60) n_e_0
}
*

**# Mortality or hospitalization analysis - Sensitivity analyses
*** sen1 - Extended the follow-up to 42 days
* Target Trial Emulation
qui forvalues k = 1/1 {
forvalues j = 0/100 {
	
	*** Prepare dataset for analysis
	capture mkdir "clinical"
	capture mkdir "clinical/sen1"

	* Setup for trial emulation
	capture mkdir "clinical/sen1/prepare"
	noi di "sen1" _col(15) "bs`j'" _col(30) "setup"
	if !fileexists("clinical/sen1/prepare/prepare sen1 bs`j'.dta") {
	use "covid paxlovid clinical main.dta", clear
	sort pseudo_key
	if `j' > 0 {
		set seed `j'
		bsample
	}

	gen date_event = min(date_death, date_admission)
	gen time_to_event = date_event - date_index
	gen treatment = inrange(time_to_treatment, 0, 1)
	
	gen censor = 1 if time_to_treatment > 5
	replace censor = 1 if date_paxlovid > date_admission
	
	gen date_censor = date_index + 5 if time_to_treatment > 5
	replace date_censor = min(date_censor, date_admission) if date_paxlovid > date_admission
	
	gen date_start_fu = date_index
	gen date_last_fu = min(date_death, mdy(02,12,2023), date_index + 42, date_molnupiravir, date_censor)
	
	gen event = inrange(date_event, date_start_fu, date_last_fu)
	gen fup_obs = min(date_event-date_index, 42) if event == 1
	replace fup_obs = min(date_last_fu-date_index, 42) if event == 0
	
	format date* %td
	* keep necessary variables
	keep pseudo_key group fup_obs event time_to_treatment time_to_event treatment censor
	gen bs = `j'
	compress
	save "clinical/sen1/prepare/prepare sen1 bs`j'.dta", replace
	}
	
	*** Cloning & censoring
	noi di "sen1" _col(15) "bs`j'" _col(30) "cloned" 
	capture mkdir "clinical/sen1/cloned"
	if !fileexists("clinical/sen1/cloned/cloned sen1 bs`j'.dta") {
	* Prepare dataset for analysis
	use "clinical/sen1/prepare/prepare sen1 bs`j'.dta", replace
	stset fup_obs, failure(event)

	* Arm A: no treatment within grace period (control: non-exposed group)
	gen outcomeA = _d // _d = event
	gen fupA = _t // _t = follow up time

	// if the patient received treatment within day 0-1:
	// 1. no event outcome, since the patient survived till censoring (treatment)
	replace outcomeA = 0 if treatment == 1 & time_to_treatment <= 1
	// 2. follow up is censored at treatment
	replace fupA = time_to_treatment if treatment == 1 & time_to_treatment <= 1

	* Arm B: treatment within grace period (treated: exposed group)
	gen outcomeB = _d 
	gen fupB = _t 

	// if the patient survived the grace period and did not receive treatment within day 0-1:
	// 1. no event outcome if the patient survived the grace period
	replace outcomeB = 0 if (treatment==0 & _t>1) | (treatment==1 & time_to_treatment >1 & time_to_treatment !=.)
	// 2. follow up is censored at end of grace period
	replace fupB = 1 if (treatment==0 & _t>1) | (treatment==1 & time_to_treatment >1 & time_to_treatment != .)

	** append clones 
	preserve
		drop outcomeB fupB
		rename outcomeA outcome
		rename fupA fup
		gen arm = "NoTreatment"
		tempfile a
		save "`a'", replace
	restore
		drop outcomeA fupA
		rename outcomeB outcome
		rename fupB fup
		gen arm = "Treatment"	
		cap append using "`a'"

	// Weight models

	sort _all
	gen NewID = _n

	** exclude those initiated at day 0 in Late arm
	drop if time_to_treatment == 0 & arm == "NoTreatment"

	** Weight model: define survival time and event indicator	
	* Early arm
	** Case 1: they do not deviate at time of treatment, but are not at risk of deviating any more
	gen wm_fup = time_to_treatment if arm == "Treatment" & time_to_treatment<=1 & time_to_treatment!=. & treatment == 1 
	gen wm_outcome = 0 if arm == "Treatment" & time_to_treatment<=1 & time_to_treatment!=. & treatment == 1 

	** Case 2: they deviate at day 1
	replace wm_fup = 1 if arm == "Treatment" & ((treatment == 0 & fup >= 1) | (time_to_treatment>1 & treatment == 1))
	replace wm_outcome = 1 if arm == "Treatment" & ((treatment == 0 & fup >= 1) | (time_to_treatment>1 & treatment == 1)) 

	** Case 3: they do not deviate, but we need to keep their survival as observed and censor them as we do not know what happens appens afterwards
	replace wm_fup = fup if arm == "Treatment" & treatment == 0 & fup < 1
	replace wm_outcome = 0 if arm == "Treatment" & treatment == 0 & fup < 1
	** add 1 days to 0-survivors
	replace wm_fup = 1 if arm == "Treatment" & wm_fup==0 

	* Late arm
	** Case 1: they do not deviate at time of treatment, but are not at risk of deviating any more
	replace wm_fup = time_to_treatment if arm == "NoTreatment" & time_to_treatment<=1 & treatment == 1 
	replace wm_outcome = 1 if arm == "NoTreatment" & time_to_treatment<=1 & treatment == 1 

	** Case 2: they deviate at day 1
	replace wm_fup = 1 if arm == "NoTreatment" & ((treatment == 0 & fup >= 1) | (time_to_treatment>1 & treatment == 1)) 
	replace wm_outcome = 0 if arm == "NoTreatment" & ((treatment == 0 & fup >= 1) | (time_to_treatment>1 & treatment == 1)) 

	** Case 3: they do not deviate, but we need to keep their survival as observed and censor them as we do not know what happens aerwards
	replace wm_fup = fup if arm == "NoTreatment" & treatment == 0 & fup < 1
	replace wm_outcome = 0 if arm == "NoTreatment" & treatment == 0 & fup < 1

	** add 1 days to 0-survivors
	replace wm_fup = 1 if arm == "NoTreatment" & wm_fup==0
	
	*** Censoring criteria - use tab to check the number of patients in each case and their outcome/censoring assignments
	* Case 1: patients admitted and initiated treatment on day 1: outcome in early, censored in late
	*tab outcome arm if time_to_treatment == 1 & time_to_event == 1
	*tab wm_outcome arm if time_to_treatment == 1 & time_to_event == 1

	* Case 2: patients admitted on day 1 with initiation after admission: outcome in late, censored in early
	*tab outcome arm if time_to_event == 1 & time_to_treatment > time_to_event 
	*tab wm_outcome arm if time_to_event == 1 & time_to_treatment > time_to_event
	replace outcome = 0 if arm == "Treatment" & time_to_event == 1 & time_to_treatment > time_to_event
	
	* Case 3: patients admitted on days 2-5 with initiation on the same day: outcome in late, censored in early
	*tab outcome arm if inrange(time_to_treatment,2,5) & time_to_treatment == time_to_event
	*tab fup arm if inrange(time_to_treatment,2,5) & time_to_treatment == time_to_event
	
	* Case 4: patients intiated treatments after 5 days: censor them at day 5
	* check all time_to_treatment > 5 censored on day 5 (or earlier)
	*tab time_to_treatment time_to_event if inrange(time_to_treatment, 0, 5) & inrange(time_to_event, 0, 5)
	*tab time_to_treatment time_to_event if inrange(time_to_event, 0, 5)
	*tab fup arm if time_to_treatment > 5 & time_to_event > 5
	*tab time_to_treatment if time_to_event <= 5
	*tab outcome arm if time_to_treatment > 5 & time_to_event > 5
	*tab wm_outcome arm if time_to_treatment > 5 & time_to_event > 5
	replace outcome = 0 if time_to_treatment > 5 & time_to_event > 5
	replace wm_outcome = 1 if time_to_treatment > 5 & time_to_event > 5
	
	gen arm_value = 1 if arm == "Treatment"
	replace arm_value = 0 if arm == "NoTreatment"
	drop arm
	order pseudo_key Early arm_value
	compress
	save "clinical/sen1/cloned/cloned sen1 bs`j'.dta", replace
	}
	
	*** Split times
	capture mkdir "clinical/sen1/split"
	noi di "sen1" _col(15) "bs`j'" _col(30) "split"
	if !fileexists("clinical/sen1/split/split sen1 bs`j'.dta") {
	use "clinical/sen1/cloned/cloned sen1 bs`j'.dta", clear
	
	** times of event
	stset fup, fail(outcome) id(NewID)
	stsplit, at(failures)

	gen tstart = _t0
	gen tstop = _t

	** times of censoring
	gen TrialEmul_cens = 1-outcome
	stset fup, fail(TrialEmul_cens) id(NewID)
	stsplit, at(failures)

	replace tstart = _t0 if tstart<_t0 & _t0 != . & _t != .
	replace tstop = _t if tstop>_t & _t0 != . & _t != .

	order pseudo_key arm tstart tstop
	sort _all
	compress
	save "clinical/sen1/split/split sen1 bs`j'.dta", replace
	}
	
	*** IPCW - Early arm
	capture mkdir "clinical/sen1/Treatment"
	noi di "sen1" _col(15) "bs`j'" _col(30) "Treatment" 
	if !fileexists("clinical/sen1/Treatment/Treatment sen1 bs`j'.dta") {
	use "clinical/sen1/split/split sen1 bs`j'.dta", clear
	keep if arm_value == 1
	merge m:1 pseudo_key using "covid paxlovid characteristics.dta", keep(3) nogen
	replace outcome = . if TrialEmul_cens == .

	* adapt records to the long format
	sort NewID tstart
	bysort NewID (tstart): replace wm_outcome = 0 if _n!=_N
	
	* Weight model:
	stset tstop, origin(time tstart) failure(wm_outcome) id(NewID)
	stcox age gender charlson_index steroid_covid_bl immuno_hist covid_hist symptomatic_bl healthcare region_num_2 region_num_3 region_num_4 vaccine_status_2 vaccine_status_3 month_covid_cate_2 month_covid_cate_3 detection_cate, efron
	if _rc == 0 {
		predict ch2_2, basech
		gen s0_2=exp(-ch2_2)
		predict xbCox2_2, xb
		gen weight = 1/s0_2^(exp(xbCox2_2)) 
		drop _st-xbCox2_2 
		gen invalid = 0
	}
	else {
		gen invalid = 1
		gen weight = .
		noi di "invalid"
	}
	keep pseudo_key arm_value tstart tstop treatment bs outcome wm_outcome NewID weight invalid
	compress
	save "clinical/sen1/Treatment/Treatment sen1 bs`j'.dta", replace
	}
	
	*** IPCW - Late arm
	capture mkdir "clinical/sen1/NoTreatment"
	noi di "sen1" _col(15) "bs`j'" _col(30) "NoTreatment" 
	if !fileexists("clinical/sen1/NoTreatment/Notreatment sen1 bs`j'.dta") {
	use "clinical/sen1/split/split sen1 bs`j'.dta", clear
	keep if arm_value == 0
	merge m:1 pseudo_key using "covid paxlovid characteristics.dta", keep(3) nogen
	replace outcome = . if TrialEmul_cens == .

	* adapt records to the long format
	sort NewID tstart
	bysort NewID: replace wm_outcome = 0 if _n!=_N

	* Weight model:
	stset tstop, origin(time tstart) failure(wm_outcome) id(NewID)
	stcox age gender charlson_index steroid_covid_bl immuno_hist covid_hist symptomatic_bl healthcare region_num_2 region_num_3 region_num_4 vaccine_status_2 vaccine_status_3 month_covid_cate_2 month_covid_cate_3 detection_cate, efron
	if _rc == 0 {
		predict ch2_2, basech
		gen s0_2=exp(-ch2_2)
		predict xbCox2_2, xb
		gen weight = 1/s0_2^(exp(xbCox2_2)) 
		drop _st-xbCox2_2 
		gen invalid = 0
	}
	else {
		gen invalid = 1
		gen weight = .
		noi di "invalid"
	}
	keep pseudo_key arm_value tstart tstop treatment bs outcome wm_outcome NewID weight invalid
	compress
	save "clinical/sen1/NoTreatment/Notreatment sen1 bs`j'.dta", replace
	}
	
	*** Combine & Generate weights
	capture mkdir "clinical/sen1/weight"
	noi di "sen1" _col(15) "bs`j'" _col(30) "Combine & weights"
	if !fileexists("clinical/sen1/weight/weight sen1 bs`j'.dta") {
	use "clinical/sen1/Treatment/Treatment sen1 bs`j'.dta", clear
	append using "clinical/sen1/NoTreatment/Notreatment sen1 bs`j'.dta"

	// create a new ID variable for each clone in each arm
	tostring NewID, gen(NewID_str)
	replace NewID_str = "00000" + NewID_str if length(NewID_str)== 1
	replace NewID_str = "0000" + NewID_str if length(NewID_str)== 2
	replace NewID_str = "000" + NewID_str if length(NewID_str)== 3
	replace NewID_str = "00" + NewID_str if length(NewID_str)== 4
	replace NewID_str = "0" + NewID_str if length(NewID_str)== 5
	gen Anal_ID = "1" + NewID_str if arm_value == 1
	replace Anal_ID = "2" + NewID_str if arm_value == 0

	replace weight = 1 if wm_outcome == 1 & tstop == 1
	sum weight
	local max_weight = r(max)

	keep pseudo_key arm_value tstart tstop treatment bs outcome wm_outcome Anal_ID weight invalid
	destring Anal_ID, replace

	* generate denominator
	rename weight _ipcw
	gen prob_denom = 1/_ipcw
	* generate numerator
	bysort arm_value tstart : gen N_group = _N
	bysort arm_value tstart : egen N_censored = sum(wm_outcome)
	gen prob_uncensored = (N_group-N_censored)/N_group
	gen prob_num = prob_uncensored
	sort arm_value pseudo_key tstart
	by arm_value pseudo_key : replace prob_num=prob_num*prob_num[_n-1] if _n!=1
	* generate stabilized IPCW
	gen _ipcw_stab = prob_num / prob_denom
	* generate combined weight
	merge m:1 pseudo_key using "covid paxlovid clinical IPTW subgp_1.dta", keep(3) keepusing(_iptw _iptw_stab)
	gen weight = _iptw * _ipcw
	gen weight_stab = _iptw_stab * _ipcw_stab
	compress
	save "clinical/sen1/weight/weight sen1 bs`j'.dta", replace
	}
	
	*** Generate KM estimate using weight_stab
	capture mkdir "clinical/sen1/KM weight_stab"
	noi di "sen1" _col(15) "bs`j'" _col(30) "KM weight_stab"
	if !fileexists("clinical/sen1/KM weight_stab/KM weight_stab sen1 bs`j'.dta") {
	use "clinical/sen1/weight/weight sen1 bs`j'.dta", clear
	count if invalid == 1
	if r(N) == 0 {
		stset tstop [pweight = weight_stab], enter(time tstart) failure(outcome)
		sts generate KM_s_w = s if arm_value == 1
		sts generate KM_ns_w = s if arm_value == 0
	}
	else {
		gen KM_s_w = .
		gen KM_ns_w = .
	}
	collapse (firstnm) KM_s_w KM_ns_w, by(tstop bs invalid)
	save "clinical/sen1/KM weight_stab/KM weight_stab sen1 bs`j'.dta", replace
	}
}
}
*
*** sen2 - Days 0-2 vs Days >2
* IPTW - clinical - Days 0-2 vs Days >2
qui forvalues k = 1/1 {
	if `k' == 1 {
		noi di "sen2" _col(15) "gp" _col(30) "mean" _col(45) "sd" _col(60) "min" _col(75) "p25" _col(90) "p50" _col(105) "p90" _col(120) "max"
	}
	
	use "covid paxlovid clinical main.dta", clear
	merge m:1 pseudo_key using "covid paxlovid characteristics.dta", keep(3) nogen

	misstable sum group age gender charlson_index steroid_covid_bl immuno_hist covid_hist symptomatic_bl healthcare region_num_2 region_num_3 region_num_4 vaccine_status_2 vaccine_status_3 month_covid_cate_2 month_covid_cate_3 detection_cate

	* paxlovid
	replace group = inrange(time_to_treatment, 0, 2)
	logit group age gender charlson_index steroid_covid_bl immuno_hist covid_hist symptomatic_bl healthcare region_num_2 region_num_3 region_num_4 vaccine_status_2 vaccine_status_3 month_covid_cate_2 month_covid_cate_3 detection_cate

	predict pr
	gen pscore = pr*group + (1-pr)*(1-group)

	* generate denominator
	count
	scalar N_all = r(N)
	bysort group : gen N_group = _N
	gen p_i = N_group / N_all
	tab p_i group

	* generate IPTW
	gen _iptw = 1 / pscore
	sum _iptw if group == 1, d
	sum _iptw if group == 0, d

	* generate stabilized IPTW
	gen _iptw_stab = p_i / pscore
	sum _iptw_stab if group == 1, d
	noi di "sen`k'" _col(15) "gp=1" _col(30) r(mean) _col(45) r(sd) _col(60) r(min) _col(75) r(p25) _col(90) r(p50) _col(105) r(p90) _col(120) r(max)
	sum _iptw_stab if group == 0, d
	noi di "sen`k'" _col(15) "gp=0" _col(30) r(mean) _col(45) r(sd) _col(60) r(min) _col(75) r(p25) _col(90) r(p50) _col(105) r(p90) _col(120) r(max)
	sort pseudo_key
	compress
	save "covid paxlovid clinical IPTW sen2.dta", replace
}
*	
* Target Trial Emulation
qui forvalues k = 1/1 {
forvalues j = 0/100 {
	
	*** Prepare dataset for analysis
	capture mkdir "clinical"
	capture mkdir "clinical/sen2"

	* Setup for trial emulation
	capture mkdir "clinical/sen2/prepare"
	noi di "sen2" _col(15) "bs`j'" _col(30) "setup"
	if !fileexists("clinical/sen2/prepare/prepare sen2 bs`j'.dta") {
	use "covid paxlovid clinical main.dta", clear
	sort pseudo_key
	if `j' > 0 {
		set seed `j'
		bsample
	}

	replace group = inrange(time_to_treatment, 0, 2)
	gen date_event = min(date_death, date_admission)
	gen time_to_event = date_event - date_index
	gen treatment = inrange(time_to_treatment, 0, 2)
	
	gen censor = 1 if time_to_treatment > 5
	replace censor = 1 if date_paxlovid > date_admission
	
	gen date_censor = date_index + 5 if time_to_treatment > 5
	replace date_censor = min(date_censor, date_admission) if date_paxlovid > date_admission
	
	gen date_start_fu = date_index
	gen date_last_fu = min(date_death, mdy(02,12,2023), date_index + 28, date_molnupiravir, date_censor)
	
	gen event = inrange(date_event, date_start_fu, date_last_fu)
	gen fup_obs = min(date_event-date_index, 28) if event == 1
	replace fup_obs = min(date_last_fu-date_index, 28) if event == 0
	
	format date* %td
	* keep necessary variables
	keep pseudo_key group fup_obs event time_to_treatment time_to_event treatment censor
	gen bs = `j'
	compress
	save "clinical/sen2/prepare/prepare sen2 bs`j'.dta", replace
	}
	
	*** Cloning & censoring
	noi di "sen2" _col(15) "bs`j'" _col(30) "cloned" 
	capture mkdir "clinical/sen2/cloned"
	if !fileexists("clinical/sen2/cloned/cloned sen2 bs`j'.dta") {
	* Prepare dataset for analysis
	use "clinical/sen2/prepare/prepare sen2 bs`j'.dta", replace
	stset fup_obs, failure(event)

	* Arm A: no treatment within grace period (control: non-exposed group)
	gen outcomeA=_d // _d=event
	gen fupA=_t // _t=follow up time

	/// if the patient received treatment within day 0-2:
	/// 1. no event outcome, since the patient survived till censoring (treatment)
	replace outcomeA=0 if treatment==1 & time_to_treatment <= 2
	/// 2. follow up is censored at treatment
	replace fupA=time_to_treatment if treatment==1 & time_to_treatment <= 2

	* Arm B: treatment within grace period (treated: exposed group)
	gen outcomeB=_d 
	gen fupB=_t 

	/// if the patient survived the grace period and did not receive treatment within day 0-2:
	/// 1. no event outcome if the patient survived the grace period
	replace outcomeB=0 if (treatment==0 & _t>2) | (treatment==1 & time_to_treatment >2 & time_to_treatment !=.)
	/// 2. follow up is censored at end of grace period
	replace fupB=1 if (treatment==0 & _t>2) | (treatment==1 & time_to_treatment >2 & time_to_treatment != .)

	** append clones 
	preserve
		drop outcomeB fupB
		rename outcomeA outcome
		rename fupA fup
		gen arm="NoTreatment"
		tempfile a
		save "`a'", replace
	restore
		drop outcomeA fupA
		rename outcomeB outcome
		rename fupB fup
		gen arm="Treatment"	
		cap append using "`a'"

	/// Weight models
	sort _all
	gen NewID=_n

	** exclude those initiated at day 0 in Late arm
	drop if time_to_treatment==0 & arm=="NoTreatment"

	** Weight model: define survival time and event indicator	
	* Early arm
	** Case 1: they do not deviate at time of treatment, but are not at risk of deviating any more
	gen wm_fup=time_to_treatment if arm=="Treatment" & time_to_treatment<=2 & time_to_treatment!=. & treatment==1 
	gen wm_outcome=0 if arm=="Treatment" & time_to_treatment<=2 & time_to_treatment!=. & treatment==1 

	** Case 2: they deviate at day 2
	replace wm_fup=2 if arm=="Treatment" & ((treatment==0 & fup>=2) | (time_to_treatment>2 & treatment==1))
	replace wm_outcome=1 if arm=="Treatment" & ((treatment==0 & fup>=2) | (time_to_treatment>2 & treatment==1)) 

	** Case 3: they do not deviate, but we need to keep their survival as observed and censor them as we do not know what happens appens afterwards
	replace wm_fup=fup if arm=="Treatment" & treatment==0 & fup<2
	replace wm_outcome=0 if arm=="Treatment" & treatment==0 & fup<2
	** add 1 days to 0-survivors
	replace wm_fup=1 if arm=="Treatment" & wm_fup==0 

	* Late arm
	** Case 1: they do not deviate at time of treatment, but are not at risk of deviating any more
	replace wm_fup=time_to_treatment if arm=="NoTreatment" & time_to_treatment<=2 & treatment==1 
	replace wm_outcome=1 if arm=="NoTreatment" & time_to_treatment<=2 & treatment==1 

	** Case 2: they deviate at day 2
	replace wm_fup=2 if arm=="NoTreatment" & ((treatment==0 & fup>=2) | (time_to_treatment>2 & treatment==1)) 
	replace wm_outcome=0 if arm=="NoTreatment" & ((treatment==0 & fup>=2) | (time_to_treatment>2 & treatment==1)) 

	** Case 3: they do not deviate, but we need to keep their survival as observed and censor them as we do not know what happens aerwards
	replace wm_fup=fup if arm=="NoTreatment" & treatment==0 & fup<2
	replace wm_outcome=0 if arm=="NoTreatment" & treatment==0 & fup<2

	** add 1 days to 0-survivors
	replace wm_fup=1 if arm=="NoTreatment" & wm_fup==0
	
	*** Censoring criteria - use tab to check the number of patients in each case and their outcome/censoring assignments
	* Case 1: patients admitted and initiated treatment on day 2: outcome in early, censored in late
	*tab outcome arm if time_to_treatment == 2 & time_to_event == 2
	*tab wm_outcome arm if time_to_treatment == 2 & time_to_event == 2

	* Case 2: patients admitted on day 2 with initiation after admission: outcome in late, censored in early
	*tab outcome arm if time_to_event == 2 & time_to_treatment > time_to_event 
	*tab wm_outcome arm if time_to_event == 2 & time_to_treatment > time_to_event
	replace outcome = 0 if arm == "Treatment" & time_to_event == 2 & time_to_treatment > time_to_event
	
	* Case 3: patients admitted on days 3-5 with initiation on the same day: outcome in late, censored in early
	*tab outcome arm if inrange(time_to_treatment,3,5) & time_to_treatment == time_to_event
	*tab fup arm if inrange(time_to_treatment,3,5) & time_to_treatment == time_to_event
	
	* Case 4: patients intiated treatments after 5 days: censor them at day 5
	* check all time_to_treatment > 5 censored on day 5 (or earlier)
	*tab time_to_treatment time_to_event if inrange(time_to_treatment, 0, 5) & inrange(time_to_event, 0, 5)
	*tab time_to_treatment time_to_event if inrange(time_to_event, 0, 5)
	*tab fup arm if time_to_treatment > 5 & time_to_event > 5
	*tab time_to_treatment if time_to_event <= 5
	*tab outcome arm if time_to_treatment > 5 & time_to_event > 5
	*tab wm_outcome arm if time_to_treatment > 5 & time_to_event > 5
	replace outcome = 0 if time_to_treatment > 5 & time_to_event > 5
	replace wm_outcome = 1 if time_to_treatment > 5 & time_to_event > 5
	
	gen arm_value = 1 if arm == "Treatment"
	replace arm_value = 0 if arm == "NoTreatment"
	drop arm
	order pseudo_key Early arm_value
	compress
	save "clinical/sen2/cloned/cloned sen2 bs`j'.dta", replace
	}
	
	*** Split times
	capture mkdir "clinical/sen2/split"
	noi di "sen2" _col(15) "bs`j'" _col(30) "split"
	if !fileexists("clinical/sen2/split/split sen2 bs`j'.dta") {
	use "clinical/sen2/cloned/cloned sen2 bs`j'.dta", clear
	
	** times of event
	stset fup, fail(outcome) id(NewID)
	stsplit, at(failures)

	gen tstart = _t0
	gen tstop = _t

	** times of censoring
	gen TrialEmul_cens = 1-outcome
	stset fup, fail(TrialEmul_cens) id(NewID)
	stsplit, at(failures)

	replace tstart = _t0 if tstart<_t0 & _t0 != . & _t != .
	replace tstop = _t if tstop>_t & _t0 != . & _t != .

	order pseudo_key arm tstart tstop
	sort _all
	compress
	save "clinical/sen2/split/split sen2 bs`j'.dta", replace
	}
	
	*** IPCW - Early arm
	capture mkdir "clinical/sen2/Treatment"
	noi di "sen2" _col(15) "bs`j'" _col(30) "Treatment" 
	if !fileexists("clinical/sen2/Treatment/Treatment sen2 bs`j'.dta") {
	use "clinical/sen2/split/split sen2 bs`j'.dta", clear
	keep if arm_value == 1
	merge m:1 pseudo_key using "covid paxlovid characteristics.dta", keep(3) nogen
	replace outcome = . if TrialEmul_cens == .

	* adapt records to the long format
	sort NewID tstart
	bysort NewID (tstart): replace wm_outcome = 0 if _n!=_N
	
	* Weight model:
	stset tstop, origin(time tstart) failure(wm_outcome) id(NewID)
	stcox age gender charlson_index steroid_covid_bl immuno_hist covid_hist symptomatic_bl healthcare region_num_2 region_num_3 region_num_4 vaccine_status_2 vaccine_status_3 month_covid_cate_2 month_covid_cate_3 detection_cate, efron
	if _rc == 0 {
		predict ch2_2, basech
		gen s0_2=exp(-ch2_2)
		predict xbCox2_2, xb
		gen weight = 1/s0_2^(exp(xbCox2_2)) 
		drop _st-xbCox2_2 
		gen invalid = 0
	}
	else {
		gen invalid = 1
		gen weight = .
		noi di "invalid"
	}
	keep pseudo_key arm_value tstart tstop treatment bs outcome wm_outcome NewID weight invalid
	compress
	save "clinical/sen2/Treatment/Treatment sen2 bs`j'.dta", replace
	}
	
	*** IPCW - Late arm
	capture mkdir "clinical/sen2/NoTreatment"
	noi di "sen2" _col(15) "bs`j'" _col(30) "NoTreatment" 
	if !fileexists("clinical/sen2/NoTreatment/Notreatment sen2 bs`j'.dta") {
	use "clinical/sen2/split/split sen2 bs`j'.dta", clear
	keep if arm_value == 0
	merge m:1 pseudo_key using "covid paxlovid characteristics.dta", keep(3) nogen
	replace outcome = . if TrialEmul_cens == .

	* adapt records to the long format
	sort NewID tstart
	bysort NewID: replace wm_outcome = 0 if _n!=_N

	* Weight model:
	stset tstop, origin(time tstart) failure(wm_outcome) id(NewID)
	stcox age gender charlson_index steroid_covid_bl immuno_hist covid_hist symptomatic_bl healthcare region_num_2 region_num_3 region_num_4 vaccine_status_2 vaccine_status_3 month_covid_cate_2 month_covid_cate_3 detection_cate, efron
	if _rc == 0 {
		predict ch2_2, basech
		gen s0_2=exp(-ch2_2)
		predict xbCox2_2, xb
		gen weight = 1/s0_2^(exp(xbCox2_2)) 
		drop _st-xbCox2_2 
		gen invalid = 0
	}
	else {
		gen invalid = 1
		gen weight = .
		noi di "invalid"
	}
	keep pseudo_key arm_value tstart tstop treatment bs outcome wm_outcome NewID weight invalid
	compress
	save "clinical/sen2/NoTreatment/Notreatment sen2 bs`j'.dta", replace
	}
	
	*** Combine & Generate weights
	capture mkdir "clinical/sen2/weight"
	noi di "sen2" _col(15) "bs`j'" _col(30) "Combine & weights"
	if !fileexists("clinical/sen2/weight/weight sen2 bs`j'.dta") {
	use "clinical/sen2/Treatment/Treatment sen2 bs`j'.dta", clear
	append using "clinical/sen2/NoTreatment/Notreatment sen2 bs`j'.dta"

	// create a new ID variable for each clone in each arm
	tostring NewID, gen(NewID_str)
	replace NewID_str = "00000" + NewID_str if length(NewID_str)== 1
	replace NewID_str = "0000" + NewID_str if length(NewID_str)== 2
	replace NewID_str = "000" + NewID_str if length(NewID_str)== 3
	replace NewID_str = "00" + NewID_str if length(NewID_str)== 4
	replace NewID_str = "0" + NewID_str if length(NewID_str)== 5
	gen Anal_ID = "1" + NewID_str if arm_value == 1
	replace Anal_ID = "2" + NewID_str if arm_value == 0

	replace weight = 1 if wm_outcome == 1 & tstop == 1
	sum weight
	local max_weight = r(max)

	keep pseudo_key arm_value tstart tstop treatment bs outcome wm_outcome Anal_ID weight invalid
	destring Anal_ID, replace

	* generate denominator
	rename weight _ipcw
	gen prob_denom = 1/_ipcw
	* generate numerator
	bysort arm_value tstart : gen N_group = _N
	bysort arm_value tstart : egen N_censored = sum(wm_outcome)
	gen prob_uncensored = (N_group-N_censored)/N_group
	gen prob_num = prob_uncensored
	sort arm_value pseudo_key tstart
	by arm_value pseudo_key : replace prob_num=prob_num*prob_num[_n-1] if _n!=1
	* generate stabilized IPCW
	gen _ipcw_stab = prob_num / prob_denom
	* generate combined weight
	merge m:1 pseudo_key using "covid paxlovid clinical IPTW sen2.dta", keep(3) keepusing(_iptw _iptw_stab)
	gen weight = _iptw * _ipcw
	gen weight_stab = _iptw_stab * _ipcw_stab
	compress
	save "clinical/sen2/weight/weight sen2 bs`j'.dta", replace
	}
	
	*** Generate KM estimate using weight_stab
	capture mkdir "clinical/sen2/KM weight_stab"
	noi di "sen2" _col(15) "bs`j'" _col(30) "KM weight_stab"
	if !fileexists("clinical/sen2/KM weight_stab/KM weight_stab sen2 bs`j'.dta") {
	use "clinical/sen2/weight/weight sen2 bs`j'.dta", clear
	count if invalid == 1
	if r(N) == 0 {
		stset tstop [pweight = weight_stab], enter(time tstart) failure(outcome)
		sts generate KM_s_w = s if arm_value == 1
		sts generate KM_ns_w = s if arm_value == 0
	}
	else {
		gen KM_s_w = .
		gen KM_ns_w = .
	}
	collapse (firstnm) KM_s_w KM_ns_w, by(tstop bs invalid)
	save "clinical/sen2/KM weight_stab/KM weight_stab sen2 bs`j'.dta", replace
	}
}
}
*
*** sen3 - Days 0-3 vs Days >3
* IPTW - clinical - Days 0-3 vs Days >3
qui forvalues k = 1/1 {
	if `k' == 1 {
		noi di "sen3" _col(15) "gp" _col(30) "mean" _col(45) "sd" _col(60) "min" _col(75) "p25" _col(90) "p50" _col(105) "p90" _col(120) "max"
	}
	
	use "covid paxlovid clinical main.dta", clear
	merge m:1 pseudo_key using "covid paxlovid characteristics.dta", keep(3) nogen

	misstable sum group age gender charlson_index steroid_covid_bl immuno_hist covid_hist symptomatic_bl healthcare region_num_2 region_num_3 region_num_4 vaccine_status_2 vaccine_status_3 month_covid_cate_2 month_covid_cate_3 detection_cate

	* paxlovid
	replace group = inrange(time_to_treatment, 0, 3)
	logit group age gender charlson_index steroid_covid_bl immuno_hist covid_hist symptomatic_bl healthcare region_num_2 region_num_3 region_num_4 vaccine_status_2 vaccine_status_3 month_covid_cate_2 month_covid_cate_3 detection_cate

	predict pr
	gen pscore = pr*group + (1-pr)*(1-group)

	* generate denominator
	count
	scalar N_all = r(N)
	bysort group : gen N_group = _N
	gen p_i = N_group / N_all
	tab p_i group

	* generate IPTW
	gen _iptw = 1 / pscore
	sum _iptw if group == 1, d
	sum _iptw if group == 0, d

	* generate stabilized IPTW
	gen _iptw_stab = p_i / pscore
	sum _iptw_stab if group == 1, d
	noi di "sen`k'" _col(15) "gp=1" _col(30) r(mean) _col(45) r(sd) _col(60) r(min) _col(75) r(p25) _col(90) r(p50) _col(105) r(p90) _col(120) r(max)
	sum _iptw_stab if group == 0, d
	noi di "sen`k'" _col(15) "gp=0" _col(30) r(mean) _col(45) r(sd) _col(60) r(min) _col(75) r(p25) _col(90) r(p50) _col(105) r(p90) _col(120) r(max)
	sort pseudo_key
	compress
	save "covid paxlovid clinical IPTW sen3.dta", replace
}
*
* Target Trial Emulation
qui forvalues k = 1/1 {
forvalues j = 0/100 {
	
	*** Prepare dataset for analysis
	capture mkdir "clinical"
	capture mkdir "clinical/sen3"

	* Setup for trial emulation
	capture mkdir "clinical/sen3/prepare"
	noi di "sen3" _col(15) "bs`j'" _col(30) "setup"
	if !fileexists("clinical/sen3/prepare/prepare sen3 bs`j'.dta") {
	use "covid paxlovid clinical main.dta", clear
	sort pseudo_key
	if `j' > 0 {
		set seed `j'
		bsample
	}

	replace group = inrange(time_to_treatment, 0, 3)
	gen date_event = min(date_death, date_admission)
	gen time_to_event = date_event - date_index
	gen treatment = inrange(time_to_treatment, 0, 3)
	
	gen censor = 1 if time_to_treatment > 5
	replace censor = 1 if date_paxlovid > date_admission
	
	gen date_censor = date_index + 5 if time_to_treatment > 5
	replace date_censor = min(date_censor, date_admission) if date_paxlovid > date_admission
	
	gen date_start_fu = date_index
	gen date_last_fu = min(date_death, mdy(02,12,2023), date_index + 28, date_molnupiravir, date_censor)
	
	gen event = inrange(date_event, date_start_fu, date_last_fu)
	gen fup_obs = min(date_event-date_index, 28) if event == 1
	replace fup_obs = min(date_last_fu-date_index, 28) if event == 0
	
	format date* %td
	* keep necessary variables
	keep pseudo_key group fup_obs event time_to_treatment time_to_event treatment censor
	gen bs = `j'
	compress
	save "clinical/sen3/prepare/prepare sen3 bs`j'.dta", replace
	}
	
	*** Cloning & censoring
	noi di "sen3" _col(15) "bs`j'" _col(30) "cloned" 
	capture mkdir "clinical/sen3/cloned"
	if !fileexists("clinical/sen3/cloned/cloned sen3 bs`j'.dta") {
	* Prepare dataset for analysis
	use "clinical/sen3/prepare/prepare sen3 bs`j'.dta", replace
	stset fup_obs, failure(event)

	* Arm A: no treatment within grace period (control: non-exposed group)
	gen outcomeA=_d // _d=event
	gen fupA=_t // _t=follow up time

	/// if the patient received treatment within day 0-1:
	/// 1. no event outcome, since the patient survived till censoring (treatment)
	replace outcomeA=0 if treatment==1 & time_to_treatment <= 3
	/// 2. follow up is censored at treatment
	replace fupA=time_to_treatment if treatment==1 & time_to_treatment <= 3

	* Arm B: treatment within grace period (treated: exposed group)
	gen outcomeB=_d 
	gen fupB=_t 

	/// if the patient survived the grace period and did not receive treatment within day 0-1:
	/// 1. no event outcome if the patient survived the grace period
	replace outcomeB=0 if (treatment==0 & _t>3) | (treatment==1 & time_to_treatment >3 & time_to_treatment !=.)
	/// 2. follow up is censored at end of grace period
	replace fupB=1 if (treatment==0 & _t>3) | (treatment==1 & time_to_treatment >3 & time_to_treatment != .)

	** append clones 
	preserve
		drop outcomeB fupB
		rename outcomeA outcome
		rename fupA fup
		gen arm="NoTreatment"
		tempfile a
		save "`a'", replace
	restore
		drop outcomeA fupA
		rename outcomeB outcome
		rename fupB fup
		gen arm="Treatment"	
		cap append using "`a'"

	/// Weight models
	sort _all
	gen NewID=_n
	
	** exclude those initiated at day 0 in Late arm
	drop if time_to_treatment==0 & arm=="NoTreatment"

	** Weight model: define survival time and event indicator	
	* Early arm
	** Case 1: they do not deviate at time of treatment, but are not at risk of deviating any more
	gen wm_fup=time_to_treatment if arm=="Treatment" & time_to_treatment<=3 & time_to_treatment!=. & treatment==1 
	gen wm_outcome=0 if arm=="Treatment" & time_to_treatment<=3 & time_to_treatment!=. & treatment==1 

	** Case 2: they deviate at day 1
	replace wm_fup=3 if arm=="Treatment" & ((treatment==0 & fup>=3) | (time_to_treatment>3 & treatment==1))
	replace wm_outcome=1 if arm=="Treatment" & ((treatment==0 & fup>=3) | (time_to_treatment>3 & treatment==1)) 

	** Case 3: they do not deviate, but we need to keep their survival as observed and censor them as we do not know what happens appens afterwards
	replace wm_fup=fup if arm=="Treatment" & treatment==0 & fup<3
	replace wm_outcome=0 if arm=="Treatment" & treatment==0 & fup<3
	** add 1 days to 0-survivors
	replace wm_fup=1 if arm=="Treatment" & wm_fup==0 

	* Late arm
	** Case 1: they do not deviate at time of treatment, but are not at risk of deviating any more
	replace wm_fup=time_to_treatment if arm=="NoTreatment" & time_to_treatment<=3 & treatment==1 
	replace wm_outcome=1 if arm=="NoTreatment" & time_to_treatment<=3 & treatment==1 

	** Case 2: they deviate at day 1
	replace wm_fup=3 if arm=="NoTreatment" & ((treatment==0 & fup>=3) | (time_to_treatment>3 & treatment==1)) 
	replace wm_outcome=0 if arm=="NoTreatment" & ((treatment==0 & fup>=3) | (time_to_treatment>3 & treatment==1)) 

	** Case 3: they do not deviate, but we need to keep their survival as observed and censor them as we do not know what happens aerwards
	replace wm_fup=fup if arm=="NoTreatment" & treatment==0 & fup<3
	replace wm_outcome=0 if arm=="NoTreatment" & treatment==0 & fup<3

	** add 1 days to 0-survivors
	replace wm_fup=1 if arm=="NoTreatment" & wm_fup==0
	
	*** Censoring criteria - use tab to check the number of patients in each case and their outcome/censoring assignments
	* Case 1: patients admitted and initiated treatment on day 1: outcome in early, censored in late
	*tab outcome arm if time_to_treatment == 3 & time_to_event == 3
	*tab wm_outcome arm if time_to_treatment == 3 & time_to_event == 3

	* Case 2: patients admitted on day 3 with initiation after admission: outcome in late, censored in early
	*tab outcome arm if time_to_event == 3 & time_to_treatment > time_to_event 
	*tab wm_outcome arm if time_to_event == 3 & time_to_treatment > time_to_event
	replace outcome = 0 if arm == "Treatment" & time_to_event == 3 & time_to_treatment > time_to_event
	
	* Case 3: patients admitted on days 4-5 with initiation on the same day: outcome in late, censored in early
	*tab outcome arm if inrange(time_to_treatment,4,5) & time_to_treatment == time_to_event
	*tab fup arm if inrange(time_to_treatment,4,5) & time_to_treatment == time_to_event
	
	* Case 4: patients intiated treatments after 5 days: censor them at day 5
	* check all time_to_treatment > 5 censored on day 5 (or earlier)
	*tab time_to_treatment time_to_event if inrange(time_to_treatment, 0, 5) & inrange(time_to_event, 0, 5)
	*tab time_to_treatment time_to_event if inrange(time_to_event, 0, 5)
	*tab fup arm if time_to_treatment > 5 & time_to_event > 5
	*tab time_to_treatment if time_to_event <= 5
	*tab outcome arm if time_to_treatment > 5 & time_to_event > 5
	*tab wm_outcome arm if time_to_treatment > 5 & time_to_event > 5
	replace outcome = 0 if time_to_treatment > 5 & time_to_event > 5
	replace wm_outcome = 1 if time_to_treatment > 5 & time_to_event > 5
	
	gen arm_value = 1 if arm == "Treatment"
	replace arm_value = 0 if arm == "NoTreatment"
	drop arm
	order pseudo_key Early arm_value
	compress
	save "clinical/sen3/cloned/cloned sen3 bs`j'.dta", replace
	}
	
	*** Split times
	capture mkdir "clinical/sen3/split"
	noi di "sen3" _col(15) "bs`j'" _col(30) "split"
	if !fileexists("clinical/sen3/split/split sen3 bs`j'.dta") {
	use "clinical/sen3/cloned/cloned sen3 bs`j'.dta", clear
	
	** times of event
	stset fup, fail(outcome) id(NewID)
	stsplit, at(failures)

	gen tstart = _t0
	gen tstop = _t

	** times of censoring
	gen TrialEmul_cens = 1-outcome
	stset fup, fail(TrialEmul_cens) id(NewID)
	stsplit, at(failures)

	replace tstart = _t0 if tstart<_t0 & _t0 != . & _t != .
	replace tstop = _t if tstop>_t & _t0 != . & _t != .

	order pseudo_key arm tstart tstop
	sort _all
	compress
	save "clinical/sen3/split/split sen3 bs`j'.dta", replace
	}
	
	*** IPCW - Early arm
	capture mkdir "clinical/sen3/Treatment"
	noi di "sen3" _col(15) "bs`j'" _col(30) "Treatment" 
	if !fileexists("clinical/sen3/Treatment/Treatment sen3 bs`j'.dta") {
	use "clinical/sen3/split/split sen3 bs`j'.dta", clear
	keep if arm_value == 1
	merge m:1 pseudo_key using "covid paxlovid characteristics.dta", keep(3) nogen
	replace outcome = . if TrialEmul_cens == .

	* adapt records to the long format
	sort NewID tstart
	bysort NewID (tstart): replace wm_outcome = 0 if _n!=_N
	
	* Weight model:
	stset tstop, origin(time tstart) failure(wm_outcome) id(NewID)
	stcox age gender charlson_index steroid_covid_bl immuno_hist covid_hist symptomatic_bl healthcare region_num_2 region_num_3 region_num_4 vaccine_status_2 vaccine_status_3 month_covid_cate_2 month_covid_cate_3 detection_cate, efron
	if _rc == 0 {
		predict ch2_2, basech
		gen s0_2=exp(-ch2_2)
		predict xbCox2_2, xb
		gen weight = 1/s0_2^(exp(xbCox2_2)) 
		drop _st-xbCox2_2 
		gen invalid = 0
	}
	else {
		gen invalid = 1
		gen weight = .
		noi di "invalid"
	}
	keep pseudo_key arm_value tstart tstop treatment bs outcome wm_outcome NewID weight invalid
	compress
	save "clinical/sen3/Treatment/Treatment sen3 bs`j'.dta", replace
	}
	
	*** IPCW - Late arm
	capture mkdir "clinical/sen3/NoTreatment"
	noi di "sen3" _col(15) "bs`j'" _col(30) "NoTreatment" 
	if !fileexists("clinical/sen3/NoTreatment/Notreatment sen3 bs`j'.dta") {
	use "clinical/sen3/split/split sen3 bs`j'.dta", clear
	keep if arm_value == 0
	merge m:1 pseudo_key using "covid paxlovid characteristics.dta", keep(3) nogen
	replace outcome = . if TrialEmul_cens == .

	* adapt records to the long format
	sort NewID tstart
	bysort NewID: replace wm_outcome = 0 if _n!=_N

	* Weight model:
	stset tstop, origin(time tstart) failure(wm_outcome) id(NewID)
	stcox age gender charlson_index steroid_covid_bl immuno_hist covid_hist symptomatic_bl healthcare region_num_2 region_num_3 region_num_4 vaccine_status_2 vaccine_status_3 month_covid_cate_2 month_covid_cate_3 detection_cate, efron
	if _rc == 0 {
		predict ch2_2, basech
		gen s0_2=exp(-ch2_2)
		predict xbCox2_2, xb
		gen weight = 1/s0_2^(exp(xbCox2_2)) 
		drop _st-xbCox2_2 
		gen invalid = 0
	}
	else {
		gen invalid = 1
		gen weight = .
		noi di "invalid"
	}
	keep pseudo_key arm_value tstart tstop treatment bs outcome wm_outcome NewID weight invalid
	compress
	save "clinical/sen3/NoTreatment/Notreatment sen3 bs`j'.dta", replace
	}
	
	*** Combine & Generate weights
	capture mkdir "clinical/sen3/weight"
	noi di "sen3" _col(15) "bs`j'" _col(30) "Combine & weights"
	if !fileexists("clinical/sen3/weight/weight sen3 bs`j'.dta") {
	use "clinical/sen3/Treatment/Treatment sen3 bs`j'.dta", clear
	append using "clinical/sen3/NoTreatment/Notreatment sen3 bs`j'.dta"

	// create a new ID variable for each clone in each arm
	tostring NewID, gen(NewID_str)
	replace NewID_str = "00000" + NewID_str if length(NewID_str)== 1
	replace NewID_str = "0000" + NewID_str if length(NewID_str)== 2
	replace NewID_str = "000" + NewID_str if length(NewID_str)== 3
	replace NewID_str = "00" + NewID_str if length(NewID_str)== 4
	replace NewID_str = "0" + NewID_str if length(NewID_str)== 5
	gen Anal_ID = "1" + NewID_str if arm_value == 1
	replace Anal_ID = "2" + NewID_str if arm_value == 0

	replace weight = 1 if wm_outcome == 1 & tstop == 1
	sum weight
	local max_weight = r(max)

	keep pseudo_key arm_value tstart tstop treatment bs outcome wm_outcome Anal_ID weight invalid
	destring Anal_ID, replace

	* generate denominator
	rename weight _ipcw
	gen prob_denom = 1/_ipcw
	* generate numerator
	bysort arm_value tstart : gen N_group = _N
	bysort arm_value tstart : egen N_censored = sum(wm_outcome)
	gen prob_uncensored = (N_group-N_censored)/N_group
	gen prob_num = prob_uncensored
	sort arm_value pseudo_key tstart
	by arm_value pseudo_key : replace prob_num=prob_num*prob_num[_n-1] if _n!=1
	* generate stabilized IPCW
	gen _ipcw_stab = prob_num / prob_denom
	* generate combined weight
	merge m:1 pseudo_key using "covid paxlovid clinical IPTW sen3.dta", keep(3) keepusing(_iptw _iptw_stab)
	gen weight = _iptw * _ipcw
	gen weight_stab = _iptw_stab * _ipcw_stab
	compress
	save "clinical/sen3/weight/weight sen3 bs`j'.dta", replace
	}
	
	*** Generate KM estimate using weight_stab
	capture mkdir "clinical/sen3/KM weight_stab"
	noi di "sen3" _col(15) "bs`j'" _col(30) "KM weight_stab"
	if !fileexists("clinical/sen3/KM weight_stab/KM weight_stab sen3 bs`j'.dta") {
	use "clinical/sen3/weight/weight sen3 bs`j'.dta", clear
	count if invalid == 1
	if r(N) == 0 {
		stset tstop [pweight = weight_stab], enter(time tstart) failure(outcome)
		sts generate KM_s_w = s if arm_value == 1
		sts generate KM_ns_w = s if arm_value == 0
	}
	else {
		gen KM_s_w = .
		gen KM_ns_w = .
	}
	collapse (firstnm) KM_s_w KM_ns_w, by(tstop bs invalid)
	save "clinical/sen3/KM weight_stab/KM weight_stab sen3 bs`j'.dta", replace
	}
}
}
*
* sen4 - Exclude those who initiated nirmatrelvir/ritonavir beyond 5 days after diagnosis or symptom onset
* IPTW - clinical - Exclude those who initiated nirmatrelvir/ritonavir beyond 5 days after diagnosis or symptom onset
qui forvalues k = 1/1 {
	if `k' == 1 {
		noi di "subgp'" _col(15) "gp" _col(30) "mean" _col(45) "sd" _col(60) "min" _col(75) "p25" _col(90) "p50" _col(105) "p90" _col(120) "max"
	}
	
	use "covid paxlovid clinical main.dta", clear
	merge m:1 pseudo_key using "covid paxlovid characteristics.dta", keep(3) nogen

	misstable sum group age gender charlson_index steroid_covid_bl immuno_hist covid_hist symptomatic_bl healthcare region_num_2 region_num_3 region_num_4 vaccine_status_2 vaccine_status_3 month_covid_cate_2 month_covid_cate_3 detection_cate

	* paxlovid
	drop if time_to_treatment > 5
	logit group age gender charlson_index steroid_covid_bl immuno_hist covid_hist symptomatic_bl healthcare region_num_2 region_num_3 region_num_4 vaccine_status_2 vaccine_status_3 month_covid_cate_2 month_covid_cate_3 detection_cate

	predict pr
	gen pscore = pr*group + (1-pr)*(1-group)

	* generate denominator
	count
	scalar N_all = r(N)
	bysort group : gen N_group = _N
	gen p_i = N_group / N_all
	tab p_i group

	* generate IPTW
	gen _iptw = 1 / pscore
	sum _iptw if group == 1, d
	sum _iptw if group == 0, d

	* generate stabilized IPTW
	gen _iptw_stab = p_i / pscore
	sum _iptw_stab if group == 1, d
	noi di "sen`k'" _col(15) "gp=1" _col(30) r(mean) _col(45) r(sd) _col(60) r(min) _col(75) r(p25) _col(90) r(p50) _col(105) r(p90) _col(120) r(max)
	sum _iptw_stab if group == 0, d
	noi di "sen`k'" _col(15) "gp=0" _col(30) r(mean) _col(45) r(sd) _col(60) r(min) _col(75) r(p25) _col(90) r(p50) _col(105) r(p90) _col(120) r(max)
	sort pseudo_key
	compress
	save "covid paxlovid clinical IPTW sen4.dta", replace
}
*
* Target Trial Emulation
qui forvalues k = 1/1 {
forvalues j = 0/100 {
	
	*** Prepare dataset for analysis
	capture mkdir "clinical"
	capture mkdir "clinical/sen4"
	
	* Setup for trial emulation
	capture mkdir "clinical/sen4/prepare"
	noi di "sen4" _col(15) "bs`j'" _col(30) "setup"
	if !fileexists("clinical/sen4/prepare/prepare sen4 bs`j'.dta") {
	use "covid paxlovid clinical main.dta", clear
	sort pseudo_key
	if `j' > 0 {
		set seed `j'
		bsample
	}

	drop if time_to_treatment > 5
	gen date_event = min(date_death, date_admission)
	gen time_to_event = date_event - date_index
	gen treatment = inrange(time_to_treatment, 0, 1)
	
	gen censor = 1 if time_to_treatment > 5
	replace censor = 1 if date_paxlovid > date_admission
	
	gen date_censor = date_index + 5 if time_to_treatment > 5
	replace date_censor = min(date_censor, date_admission) if date_paxlovid > date_admission
	
	gen date_start_fu = date_index
	gen date_last_fu = min(date_death, mdy(02,12,2023), date_index + 28, date_molnupiravir, date_censor)
	
	gen event = inrange(date_event, date_start_fu, date_last_fu)
	gen fup_obs = min(date_event-date_index, 28) if event == 1
	replace fup_obs = min(date_last_fu-date_index, 28) if event == 0
	
	format date* %td
	* keep necessary variables
	keep pseudo_key group fup_obs event time_to_treatment time_to_event treatment censor
	gen bs = `j'
	compress
	save "clinical/sen4/prepare/prepare sen4 bs`j'.dta", replace
	}
	
	*** Cloning & censoring
	noi di "sen4" _col(15) "bs`j'" _col(30) "cloned" 
	capture mkdir "clinical/sen4/cloned"
	if !fileexists("clinical/sen4/cloned/cloned sen4 bs`j'.dta") {
	* Prepare dataset for analysis
	use "clinical/sen4/prepare/prepare sen4 bs`j'.dta", replace
	stset fup_obs, failure(event)

	* Arm A: no treatment within grace period (control: non-exposed group)
	gen outcomeA = _d // _d = event
	gen fupA = _t // _t = follow up time

	/// if the patient received treatment within day 0-1:
	/// 1. no event outcome, since the patient survived till censoring (treatment)
	replace outcomeA = 0 if treatment == 1 & time_to_treatment <= 1
	/// 2. follow up is censored at treatment
	replace fupA = time_to_treatment if treatment == 1 & time_to_treatment <= 1

	* Arm B: treatment within grace period (treated: exposed group)
	gen outcomeB = _d 
	gen fupB = _t 

	/// if the patient survived the grace period and did not receive treatment within day 0-1:
	/// 1. no event outcome if the patient survived the grace period
	replace outcomeB = 0 if (treatment==0 & _t>1) | (treatment==1 & time_to_treatment >1 & time_to_treatment !=.)
	/// 2. follow up is censored at end of grace period
	replace fupB = 1 if (treatment==0 & _t>1) | (treatment==1 & time_to_treatment >1 & time_to_treatment != .)

	** append clones 
	preserve
		drop outcomeB fupB
		rename outcomeA outcome
		rename fupA fup
		gen arm = "NoTreatment"
		tempfile a
		save "`a'", replace
	restore
		drop outcomeA fupA
		rename outcomeB outcome
		rename fupB fup
		gen arm = "Treatment"	
		cap append using "`a'"

	/// Weight models
	sort _all
	gen NewID = _n
	
	** exclude those initiated at day 0 in Late arm
	drop if time_to_treatment == 0 & arm == "NoTreatment"

	** Weight model: define survival time and event indicator	
	* Early arm
	** Case 1: they do not deviate at time of treatment, but are not at risk of deviating any more
	gen wm_fup = time_to_treatment if arm == "Treatment" & time_to_treatment<=1 & time_to_treatment!=. & treatment == 1 
	gen wm_outcome = 0 if arm == "Treatment" & time_to_treatment<=1 & time_to_treatment!=. & treatment == 1 

	** Case 2: they deviate at day 1
	replace wm_fup = 1 if arm == "Treatment" & ((treatment == 0 & fup >= 1) | (time_to_treatment>1 & treatment == 1))
	replace wm_outcome = 1 if arm == "Treatment" & ((treatment == 0 & fup >= 1) | (time_to_treatment>1 & treatment == 1)) 

	** Case 3: they do not deviate, but we need to keep their survival as observed and censor them as we do not know what happens appens afterwards
	replace wm_fup = fup if arm == "Treatment" & treatment == 0 & fup < 1
	replace wm_outcome = 0 if arm == "Treatment" & treatment == 0 & fup < 1
	** add 1 days to 0-survivors
	replace wm_fup = 1 if arm == "Treatment" & wm_fup==0 

	* Late arm
	** Case 1: they do not deviate at time of treatment, but are not at risk of deviating any more
	replace wm_fup = time_to_treatment if arm == "NoTreatment" & time_to_treatment<=1 & treatment == 1 
	replace wm_outcome = 1 if arm == "NoTreatment" & time_to_treatment<=1 & treatment == 1 

	** Case 2: they deviate at day 1
	replace wm_fup = 1 if arm == "NoTreatment" & ((treatment == 0 & fup >= 1) | (time_to_treatment>1 & treatment == 1)) 
	replace wm_outcome = 0 if arm == "NoTreatment" & ((treatment == 0 & fup >= 1) | (time_to_treatment>1 & treatment == 1)) 

	** Case 3: they do not deviate, but we need to keep their survival as observed and censor them as we do not know what happens aerwards
	replace wm_fup = fup if arm == "NoTreatment" & treatment == 0 & fup < 1
	replace wm_outcome = 0 if arm == "NoTreatment" & treatment == 0 & fup < 1

	** add 1 days to 0-survivors
	replace wm_fup = 1 if arm == "NoTreatment" & wm_fup==0
	
	*** Censoring criteria - use tab to check the number of patients in each case and their outcome/censoring assignments
	* Case 1: patients admitted and initiated treatment on day 1: outcome in early, censored in late
	*tab outcome arm if time_to_treatment == 1 & time_to_event == 1
	*tab wm_outcome arm if time_to_treatment == 1 & time_to_event == 1

	* Case 2: patients admitted on day 1 with initiation after admission: outcome in late, censored in early
	*tab outcome arm if time_to_event == 1 & time_to_treatment > time_to_event 
	*tab wm_outcome arm if time_to_event == 1 & time_to_treatment > time_to_event
	replace outcome = 0 if arm == "Treatment" & time_to_event == 1 & time_to_treatment > time_to_event
	
	* Case 3: patients admitted on days 2-5 with initiation on the same day: outcome in late, censored in early
	*tab outcome arm if inrange(time_to_treatment,2,5) & time_to_treatment == time_to_event
	*tab fup arm if inrange(time_to_treatment,2,5) & time_to_treatment == time_to_event
	
	* Case 4: patients intiated treatments after 5 days: censor them at day 5
	* check all time_to_treatment > 5 censored on day 5 (or earlier)
	*tab time_to_treatment time_to_event if inrange(time_to_treatment, 0, 5) & inrange(time_to_event, 0, 5)
	*tab time_to_treatment time_to_event if inrange(time_to_event, 0, 5)
	*tab fup arm if time_to_treatment > 5 & time_to_event > 5
	*tab time_to_treatment if time_to_event <= 5
	*tab outcome arm if time_to_treatment > 5 & time_to_event > 5
	*tab wm_outcome arm if time_to_treatment > 5 & time_to_event > 5
	replace outcome = 0 if time_to_treatment > 5 & time_to_event > 5
	replace wm_outcome = 1 if time_to_treatment > 5 & time_to_event > 5
	
	gen arm_value = 1 if arm == "Treatment"
	replace arm_value = 0 if arm == "NoTreatment"
	drop arm
	order pseudo_key Early arm_value
	compress
	save "clinical/sen4/cloned/cloned sen4 bs`j'.dta", replace
	}
	
	*** Split times
	capture mkdir "clinical/sen4/split"
	noi di "sen4" _col(15) "bs`j'" _col(30) "split"
	if !fileexists("clinical/sen4/split/split sen4 bs`j'.dta") {
	use "clinical/sen4/cloned/cloned sen4 bs`j'.dta", clear
	
	** times of event
	stset fup, fail(outcome) id(NewID)
	stsplit, at(failures)

	gen tstart = _t0
	gen tstop = _t

	** times of censoring
	gen TrialEmul_cens = 1-outcome
	stset fup, fail(TrialEmul_cens) id(NewID)
	stsplit, at(failures)

	replace tstart = _t0 if tstart<_t0 & _t0 != . & _t != .
	replace tstop = _t if tstop>_t & _t0 != . & _t != .

	order pseudo_key arm tstart tstop
	sort _all
	compress
	save "clinical/sen4/split/split sen4 bs`j'.dta", replace
	}
	
	*** IPCW - Early arm
	capture mkdir "clinical/sen4/Treatment"
	noi di "sen4" _col(15) "bs`j'" _col(30) "Treatment" 
	if !fileexists("clinical/sen4/Treatment/Treatment sen4 bs`j'.dta") {
	use "clinical/sen4/split/split sen4 bs`j'.dta", clear
	keep if arm_value == 1
	merge m:1 pseudo_key using "covid paxlovid characteristics.dta", keep(3) nogen
	replace outcome = . if TrialEmul_cens == .

	* adapt records to the long format
	sort NewID tstart
	bysort NewID (tstart): replace wm_outcome = 0 if _n!=_N
	
	* Weight model:
	stset tstop, origin(time tstart) failure(wm_outcome) id(NewID)
	stcox age gender charlson_index steroid_covid_bl immuno_hist covid_hist symptomatic_bl healthcare region_num_2 region_num_3 region_num_4 vaccine_status_2 vaccine_status_3 month_covid_cate_2 month_covid_cate_3 detection_cate, efron
	if _rc == 0 {
		predict ch2_2, basech
		gen s0_2=exp(-ch2_2)
		predict xbCox2_2, xb
		gen weight = 1/s0_2^(exp(xbCox2_2)) 
		drop _st-xbCox2_2 
		gen invalid = 0
	}
	else {
		gen invalid = 1
		gen weight = .
		noi di "invalid"
	}
	keep pseudo_key arm_value tstart tstop treatment bs outcome wm_outcome NewID weight invalid
	compress
	save "clinical/sen4/Treatment/Treatment sen4 bs`j'.dta", replace
	}
	
	*** IPCW - Late arm
	capture mkdir "clinical/sen4/NoTreatment"
	noi di "sen4" _col(15) "bs`j'" _col(30) "NoTreatment" 
	if !fileexists("clinical/sen4/NoTreatment/Notreatment sen4 bs`j'.dta") {
	use "clinical/sen4/split/split sen4 bs`j'.dta", clear
	keep if arm_value == 0
	merge m:1 pseudo_key using "covid paxlovid characteristics.dta", keep(3) nogen
	replace outcome = . if TrialEmul_cens == .

	* adapt records to the long format
	sort NewID tstart
	bysort NewID: replace wm_outcome = 0 if _n!=_N

	* Weight model:
	stset tstop, origin(time tstart) failure(wm_outcome) id(NewID)
	stcox age gender charlson_index steroid_covid_bl immuno_hist covid_hist symptomatic_bl healthcare region_num_2 region_num_3 region_num_4 vaccine_status_2 vaccine_status_3 month_covid_cate_2 month_covid_cate_3 detection_cate, efron
	if _rc == 0 {
		predict ch2_2, basech
		gen s0_2=exp(-ch2_2)
		predict xbCox2_2, xb
		gen weight = 1/s0_2^(exp(xbCox2_2)) 
		drop _st-xbCox2_2 
		gen invalid = 0
	}
	else {
		gen invalid = 1
		gen weight = .
		noi di "invalid"
	}
	keep pseudo_key arm_value tstart tstop treatment bs outcome wm_outcome NewID weight invalid
	compress
	save "clinical/sen4/NoTreatment/Notreatment sen4 bs`j'.dta", replace
	}
	
	*** Combine & Generate weights
	capture mkdir "clinical/sen4/weight"
	noi di "sen4" _col(15) "bs`j'" _col(30) "Combine & weights"
	if !fileexists("clinical/sen4/weight/weight sen4 bs`j'.dta") {
	use "clinical/sen4/Treatment/Treatment sen4 bs`j'.dta", clear
	append using "clinical/sen4/NoTreatment/Notreatment sen4 bs`j'.dta"

	// create a new ID variable for each clone in each arm
	tostring NewID, gen(NewID_str)
	replace NewID_str = "00000" + NewID_str if length(NewID_str)== 1
	replace NewID_str = "0000" + NewID_str if length(NewID_str)== 2
	replace NewID_str = "000" + NewID_str if length(NewID_str)== 3
	replace NewID_str = "00" + NewID_str if length(NewID_str)== 4
	replace NewID_str = "0" + NewID_str if length(NewID_str)== 5
	gen Anal_ID = "1" + NewID_str if arm_value == 1
	replace Anal_ID = "2" + NewID_str if arm_value == 0

	replace weight = 1 if wm_outcome == 1 & tstop == 1
	sum weight
	local max_weight = r(max)

	keep pseudo_key arm_value tstart tstop treatment bs outcome wm_outcome Anal_ID weight invalid
	destring Anal_ID, replace

	* generate denominator
	rename weight _ipcw
	gen prob_denom = 1/_ipcw
	* generate numerator
	bysort arm_value tstart : gen N_group = _N
	bysort arm_value tstart : egen N_censored = sum(wm_outcome)
	gen prob_uncensored = (N_group-N_censored)/N_group
	gen prob_num = prob_uncensored
	sort arm_value pseudo_key tstart
	by arm_value pseudo_key : replace prob_num=prob_num*prob_num[_n-1] if _n!=1
	* generate stabilized IPCW
	gen _ipcw_stab = prob_num / prob_denom
	* generate combined weight
	merge m:1 pseudo_key using "covid paxlovid clinical IPTW sen4.dta", keep(3) keepusing(_iptw _iptw_stab)
	gen weight = _iptw * _ipcw
	gen weight_stab = _iptw_stab * _ipcw_stab
	compress
	save "clinical/sen4/weight/weight sen4 bs`j'.dta", replace
	}
	
	*** Generate KM estimate using weight_stab
	capture mkdir "clinical/sen4/KM weight_stab"
	noi di "sen4" _col(15) "bs`j'" _col(30) "KM weight_stab"
	if !fileexists("clinical/sen4/KM weight_stab/KM weight_stab sen4 bs`j'.dta") {
	use "clinical/sen4/weight/weight sen4 bs`j'.dta", clear
	count if invalid == 1
	if r(N) == 0 {
		stset tstop [pweight = weight_stab], enter(time tstart) failure(outcome)
		sts generate KM_s_w = s if arm_value == 1
		sts generate KM_ns_w = s if arm_value == 0
	}
	else {
		gen KM_s_w = .
		gen KM_ns_w = .
	}
	collapse (firstnm) KM_s_w KM_ns_w, by(tstop bs invalid)
	save "clinical/sen4/KM weight_stab/KM weight_stab sen4 bs`j'.dta", replace
	}
}
}
*

cls
* Finalize bootstrap datasets
qui forvalues k = 1/18 {
clear
forvalues j = 0/0 {
	capture append using "clinical/sen`k'/KM weight_stab/KM weight_stab sen`k' bs`j'.dta"
}
	gen hazard_s_w = 1 - KM_s_w
	gen hazard_ns_w = 1 - KM_ns_w
	gen odds_s_w = hazard_s_w / (1-hazard_s_w)
	gen odds_ns_w = hazard_ns_w / (1-hazard_ns_w)
	gen RR_w = hazard_s_w / hazard_ns_w
	gen diff_w = hazard_s_w - hazard_ns_w
	gen OR_w = odds_s_w / odds_ns_w
	rename tstop fup
	compress
	save "clinical/sen`k'/KM weight_stab sen`k' bs_all.dta", replace
}
*
* KM estimate
qui forvalues k = 1/18 {
	use "clinical/sen`k'/KM weight_stab sen`k' bs_all.dta", clear
	
	keep if KM_s_w < . & KM_ns_w < .
	drop if invalid == 1
	bysort bs (fup) : keep if _n == _N
	keep if _n <= 1 + 100
	
	sum hazard_s_w if bs == 0, d
	scalar hazard_s_mean = r(mean)
	centile hazard_s_w if bs > 0, centile(2.5 97.5)
	scalar hazard_s_cil = r(c_1)
	scalar hazard_s_ciu = r(c_2)
	
	sum hazard_ns_w if bs == 0, d
	scalar hazard_ns_mean = r(mean)
	centile hazard_ns_w if bs > 0, centile(2.5 97.5)
	scalar hazard_ns_cil = r(c_1)
	scalar hazard_ns_ciu = r(c_2)
	
	noi di "sen`k'" _col(15) hazard_s_mean _col(30) hazard_s_cil _col(45) hazard_s_ciu _col(60) ///
	hazard_ns_mean _col(75) hazard_s_cil _col(90) hazard_s_ciu
}
*
* Absolute risk reduction
qui forvalues k = 1/18 {
	use "clinical/sen`k'/KM weight_stab sen`k' bs_all.dta", clear
	
	keep if KM_s_w < . & KM_ns_w < .
	drop if invalid == 1
	bysort bs (fup) : keep if _n == _N
	keep if _n <= 1 + 100
	
	sum diff_w if bs == 0, d
	scalar bs_mean = r(mean)
	scalar bs_p50 = r(p50)
	
	centile diff_w if bs > 0, centile(2.5 97.5)
	scalar bs_cil = r(c_1)
	scalar bs_ciu = r(c_2)
	noi di "sen`k'" _col(15) bs_mean _col(30) bs_p50 _col(45) bs_cil _col(60) bs_ciu
}
*
* Relative risk
qui forvalues k = 1/18 {
	use "clinical/sen`k'/KM weight_stab sen`k' bs_all.dta", clear
	
	keep if KM_s_w < . & KM_ns_w < .
	drop if invalid == 1
	bysort bs (fup) : keep if _n == _N
	keep if _n <= 1 + 100
	
	sum RR_w if bs == 0, d
	scalar bs_mean = r(mean)
	scalar bs_p50 = r(p50)
	centile RR_w if bs > 0, centile(2.5 97.5)
	scalar bs_cil = r(c_1)
	scalar bs_ciu = r(c_2)
	noi di "sen`k'" _col(15) bs_mean _col(30) bs_p50 _col(45) bs_cil _col(60) bs_ciu
}
*
* N
qui forvalues k = 1/18 {
	if `k' == 1 {
		noi di "subgp" _col(15) "N_1" _col(30) "N_0" _col(45) "n_e_1" _col(60) "n_e_0"
	}
	use "clinical/sen`k'/cloned/cloned sen`k' bs0.dta", clear
	*keep if (arm_value == 1 & group == 1) | (arm_value == 0 & group == 0)
	count if arm_value == 1
	scalar N_1 = r(N)
	count if arm_value == 1 & outcome == 1
	scalar n_e_1 = r(N)
	count if arm_value == 0
	scalar N_0 = r(N)
	count if arm_value == 0 & outcome == 1
	scalar n_e_0 = r(N)
	noi di "sen`k'" _col(15) N_1 _col(30) N_0 _col(45) n_e_1 _col(60) n_e_0
}
*
qui forvalues k = 1/18 {
	use "clinical/sen`k'/cloned/cloned sen`k' bs0.dta", clear
	keep if (arm_value == 1 & group == 1) | (arm_value == 0 & group == 0)
	count if group == 1
	scalar N_1 = r(N)
	count if group == 0
	scalar N_0 = r(N)
	count if group == 1 & outcome == 1
	scalar n_e_1 = r(N)
	count if group == 0 & outcome == 1
	scalar n_e_0 = r(N)
	noi di "sen`k'" _col(15) N_1 _col(30) N_0 _col(45) n_e_1 _col(60) n_e_0
}
*

**# Table 3. Association between timing of nirmatrelvir/ritonavir initiation and 28-day viral burden rebound
cls
**# Viral burden rebound analysis - Main analysis
qui forvalues k = 1/18 {
	forvalues j = 0/100 {

	capture mkdir "VBR"
	capture mkdir "VBR/subgp_`k'"
	capture mkdir "VBR/subgp_`k'/prepare"
	capture mkdir "VBR/subgp_`k'/cloned"
	capture mkdir "VBR/subgp_`k'/split"
	capture mkdir "VBR/subgp_`k'/Treatment"
	capture mkdir "VBR/subgp_`k'/NoTreatment"
	capture mkdir "VBR/subgp_`k'/weight"
	capture mkdir "VBR/subgp_`k'/KM weight_stab"
	
*** Prepare dataset for analysis
	noi di "subgp_`k'" _col(15) "bs`j'" _col(30) "setup"
	* Setup for trial emulation
	if !fileexists("VBR/subgp_`k'/prepare/prepare subgp_`k' bs`j'.dta") {
	use "covid paxlovid VBR main.dta", clear
	
	* generate subgroups
	gen subgp_1 = 1

	gen subgp_2 = gender == 1
	gen subgp_3 = gender == 0

	gen subgp_4 = month_covid_cate == 1
	gen subgp_5 = month_covid_cate == 2
	gen subgp_6 = month_covid_cate == 3

	gen subgp_7 = vaccine_status > 1
	gen subgp_8 = vaccine_status == 1

	gen subgp_9 = inrange(charlson_index, 0, 6)
	gen subgp_10 = inrange(charlson_index, 7, 12)

	gen subgp_11 = date_steroid_covid <= date_index
	gen subgp_12 = !(date_steroid_covid <= date_index)

	gen subgp_13 = immuno_hist == 1
	gen subgp_14 = immuno_hist == 0

	gen subgp_15 = date_index == date_onset
	gen subgp_16 = date_index == date_covid
	
	gen subgp_17 = inpatient_paxlovid == 1
	gen subgp_18 = inpatient_paxlovid == 0

	keep if subgp_`k' == 1
	sort pseudo_key
	if `j' > 0 {
		set seed `j'
		bsample
	}
	gen date_event = date_VBR_end
	
	* VBR
	* organize / rename / generate variables
	gen time_to_event = date_event - date_index
	gen treatment = inrange(time_to_treatment, 0, 1)
	
	gen censor = 1 if time_to_treatment > 5
	gen date_censor = date_index + 5 if time_to_treatment > 5
	
	gen date_start_fu = date_index
	gen date_last_fu = min(date_death, mdy(02,12,2023), date_start_fu + 28, date_molnupiravir, date_censor)
	
	gen event = inrange(date_event, date_start_fu, date_last_fu)
	gen fup_obs = min(date_event-date_index, 28) if event == 1
	replace fup_obs = min(date_last_fu-date_index, 28) if event == 0

	* keep necessary variables
	keep pseudo_key group fup_obs event time_to_treatment time_to_event treatment censor
	gen bs = `j'
	compress
	save "VBR/subgp_`k'/prepare/prepare subgp_`k' bs`j'.dta", replace
	}
	
	*** Cloning & censoring
	noi di "subgp_`k'" _col(15) "bs`j'" _col(30) "cloning"
	if !fileexists("VBR/subgp_`k'/cloned/cloned subgp_`k' bs`j'.dta") {
	* Prepare dataset for analysis
	use "VBR/subgp_`k'/prepare/prepare subgp_`k' bs`j'.dta", replace
	stset fup_obs, failure(event)

	* Arm A: no treatment within grace period (control: non-exposed group)
	gen outcomeA = _d // _d = event
	gen fupA = _t // _t = follow up time

	/// if the patient received treatment within day 0-1:
	/// 1. no event outcome, since the patient survived till censoring (treatment)
	replace outcomeA = 0 if treatment == 1 & time_to_treatment <= 1
	/// 2. follow up is censored at treatment
	replace fupA = time_to_treatment if treatment == 1 & time_to_treatment <= 1

	* Arm B: treatment within grace period (treated: exposed group)
	gen outcomeB = _d 
	gen fupB = _t 

	/// if the patient survived the grace period and did not receive treatment within day 0-1:
	/// 1. no event outcome if the patient survived the grace period
	replace outcomeB = 0 if (treatment==0 & _t>1) | (treatment==1 & time_to_treatment >1 & time_to_treatment !=.)
	/// 2. follow up is censored at end of grace period
	replace fupB = 1 if (treatment==0 & _t>1) | (treatment==1 & time_to_treatment >1 & time_to_treatment != .)

	** append clones 
	preserve
		drop outcomeB fupB
		rename outcomeA outcome
		rename fupA fup
		gen arm = "NoTreatment"
		tempfile a
		save "`a'", replace
	restore
		drop outcomeA fupA
		rename outcomeB outcome
		rename fupB fup
		gen arm = "Treatment"	
		cap append using "`a'"

	// Weight models
	sort _all
	gen NewID = _n

	** exclude those initiated at day 0 in Late arm
	drop if time_to_treatment == 0 & arm == "NoTreatment"

	** Weight model: define survival time and event indicator	
	* Early arm
	** Case 1: they do not deviate at time of treatment, but are not at risk of deviating any more
	gen wm_fup = time_to_treatment if arm == "Treatment" & time_to_treatment<=1 & time_to_treatment!=. & treatment == 1 
	gen wm_outcome = 0 if arm == "Treatment" & time_to_treatment<=1 & time_to_treatment!=. & treatment == 1 

	** Case 2: they deviate at day 1
	replace wm_fup = 1 if arm == "Treatment" & ((treatment == 0 & fup >= 1) | (time_to_treatment>1 & treatment == 1))
	replace wm_outcome = 1 if arm == "Treatment" & ((treatment == 0 & fup >= 1) | (time_to_treatment>1 & treatment == 1)) 

	** Case 3: they do not deviate, but we need to keep their survival as observed and censor them as we do not know what happens appens afterwards
	replace wm_fup = fup if arm == "Treatment" & treatment == 0 & fup < 1
	replace wm_outcome = 0 if arm == "Treatment" & treatment == 0 & fup < 1
	** add 1 days to 0-survivors
	replace wm_fup= 1 if arm == "Treatment" & wm_fup==0 

	* Late arm
	** Case 1: they do not deviate at time of treatment, but are not at risk of deviating any more
	replace wm_fup = time_to_treatment if arm == "NoTreatment" & time_to_treatment<=1 & treatment == 1 
	replace wm_outcome = 1 if arm == "NoTreatment" & time_to_treatment<=1 & treatment == 1 

	** Case 2: they deviate at day 1
	replace wm_fup = 1 if arm == "NoTreatment" & ((treatment == 0 & fup >= 1) | (time_to_treatment>1 & treatment == 1)) 
	replace wm_outcome = 0 if arm == "NoTreatment" & ((treatment == 0 & fup >= 1) | (time_to_treatment>1 & treatment == 1)) 

	** Case 3: they do not deviate, but we need to keep their survival as observed and censor them as we do not know what happens aerwards
	replace wm_fup = fup if arm == "NoTreatment" & treatment == 0 & fup < 1
	replace wm_outcome = 0 if arm == "NoTreatment" & treatment == 0 & fup < 1

	** add 1 days to 0-survivors
	replace wm_fup = 1 if arm == "NoTreatment" & wm_fup==0

	** censor patient at day 5 if initiated after day 5
	replace wm_outcome = 1 if time_to_treatment > 5
	
	gen arm_value = 1 if arm == "Treatment"
	replace arm_value = 0 if arm == "NoTreatment"
	drop arm
	order pseudo_key Early arm_value
	compress
	save "VBR/subgp_`k'/cloned/cloned subgp_`k' bs`j'.dta", replace
	}
	
* Split times
	noi di "subgp_`k'" _col(15) "bs`j'" _col(30) "split"
	if !fileexists("VBR/subgp_`k'/split/split subgp_`k' bs`j'.dta") {
	use "VBR/subgp_`k'/cloned/cloned subgp_`k' bs`j'.dta", clear
	
	** times of event
	stset fup, fail(outcome) id(NewID)
	stsplit, at(failures)

	gen tstart = _t0
	gen tstop = _t

	** times of censoring
	gen TrialEmul_cens = 1-outcome
	stset fup, fail(TrialEmul_cens) id(NewID)
	stsplit, at(failures)

	replace tstart = _t0 if tstart<_t0 & _t0 != . & _t != .
	replace tstop = _t if tstop>_t & _t0 != . & _t != .

	order pseudo_key arm tstart tstop
	sort _all
	compress
	save "VBR/subgp_`k'/split/split subgp_`k' bs`j'.dta", replace
	}
	
* IPCW - Early arm
	noi di "subgp_`k'" _col(15) "bs`j'" _col(30) "Treatment"
	if !fileexists("VBR/subgp_`k'/Treatment/Treatment subgp_`k' bs`j'.dta") {
	use "VBR/subgp_`k'/split/split subgp_`k' bs`j'.dta", clear
	keep if arm_value == 1
	merge m:1 pseudo_key using "covid paxlovid characteristics.dta", keep(3) nogen
	replace outcome = . if TrialEmul_cens == .

	* adapt records to the long format
	sort NewID tstart
	bysort NewID (tstart): replace wm_outcome = 0 if _n!=_N
	
	* Weight model:
	stset tstop, origin(time tstart) failure(wm_outcome) id(NewID)
	stcox age gender charlson_index steroid_covid_bl immuno_hist covid_hist symptomatic_bl healthcare region_num_2 region_num_3 region_num_4 vaccine_status_2 vaccine_status_3 month_covid_cate_2 month_covid_cate_3 detection_cate inpatient_paxlovid censor_ctv_day_5, efron
	if _rc == 0 {
		predict ch2_2, basech
		gen s0_2=exp(-ch2_2)
		predict xbCox2_2, xb
		gen weight = 1/s0_2^(exp(xbCox2_2)) 
		drop _st-xbCox2_2 
		gen invalid = 0
	}
	else {
		gen invalid = 1
		gen weight = .
		noi di "invalid"
	}
	keep pseudo_key arm_value tstart tstop treatment bs outcome wm_outcome NewID weight invalid
	compress
	save "VBR/subgp_`k'/Treatment/Treatment subgp_`k' bs`j'.dta", replace
}

* IPCW - Late arm
	noi di "subgp_`k'" _col(15) "bs`j'" _col(30) "NoTreatment"
	if !fileexists("VBR/subgp_`k'/NoTreatment/Notreatment subgp_`k' bs`j'.dta") {
	use "VBR/subgp_`k'/split/split subgp_`k' bs`j'.dta", clear
	keep if arm_value == 0
	merge m:1 pseudo_key using "covid paxlovid characteristics.dta", keep(3) nogen
	replace outcome = . if TrialEmul_cens == .

	* adapt records to the long format
	sort NewID tstart
	bysort NewID: replace wm_outcome = 0 if _n!=_N

	* Weight model:
	stset tstop, origin(time tstart) failure(wm_outcome) id(NewID)
	stcox age gender charlson_index steroid_covid_bl immuno_hist covid_hist symptomatic_bl healthcare region_num_2 region_num_3 region_num_4 vaccine_status_2 vaccine_status_3 month_covid_cate_2 month_covid_cate_3 detection_cate inpatient_paxlovid censor_ctv_day_5, efron
	if _rc == 0 {
		predict ch2_2, basech
		gen s0_2=exp(-ch2_2)
		predict xbCox2_2, xb
		gen weight = 1/s0_2^(exp(xbCox2_2)) 
		drop _st-xbCox2_2 
		gen invalid = 0
	}
	else {
		gen invalid = 1
		gen weight = .
		noi di "invalid"
	}
	keep pseudo_key arm_value tstart tstop treatment bs outcome wm_outcome NewID weight invalid
	compress
	save "VBR/subgp_`k'/NoTreatment/Notreatment subgp_`k' bs`j'.dta", replace
	}
	
* Combine & Generate weights
	noi di "subgp_`k'" _col(15) "bs`j'" _col(30) "Combine"
	if !fileexists("VBR/subgp_`k'/weight/weight subgp_`k' bs`j'.dta") {
	use "VBR/subgp_`k'/Treatment/Treatment subgp_`k' bs`j'.dta", clear
	append using "VBR/subgp_`k'/NoTreatment/Notreatment subgp_`k' bs`j'.dta"

	// create a new ID variable for each clone in each arm
	tostring NewID, gen(NewID_str)
	replace NewID_str = "00000" + NewID_str if length(NewID_str)== 1
	replace NewID_str = "0000" + NewID_str if length(NewID_str)== 2
	replace NewID_str = "000" + NewID_str if length(NewID_str)== 3
	replace NewID_str = "00" + NewID_str if length(NewID_str)== 4
	replace NewID_str = "0" + NewID_str if length(NewID_str)== 5
	gen Anal_ID = "1" + NewID_str if arm_value == 1
	replace Anal_ID = "2" + NewID_str if arm_value == 0

	replace weight = 1 if wm_outcome == 1 & tstop == 1
	sum weight
	local max_weight = r(max)
	
	keep pseudo_key arm_value tstart tstop treatment bs outcome wm_outcome Anal_ID weight invalid
	destring Anal_ID, replace

	* generate denominator
	rename weight _ipcw
	gen prob_denom = 1/_ipcw
	* generate numerator
	bysort arm_value tstart : gen N_group = _N
	bysort arm_value tstart : egen N_censored = sum(wm_outcome)
	gen prob_uncensored = (N_group-N_censored)/N_group
	gen prob_num = prob_uncensored
	sort arm_value pseudo_key tstart
	by arm_value pseudo_key : replace prob_num=prob_num*prob_num[_n-1] if _n!=1
	* generate stabilized IPCW
	gen _ipcw_stab = prob_num / prob_denom
	* generate combined weight
	merge m:1 pseudo_key using "covid paxlovid VBR IPTW subgp_`k'.dta", keep(3) keepusing(_iptw _iptw_stab)
	gen weight = _iptw * _ipcw
	gen weight_stab = _iptw_stab * _ipcw_stab
	compress
	save "VBR/subgp_`k'/weight/weight subgp_`k' bs`j'.dta", replace
	}
	
* Generate KM estimate
	noi di "subgp_`k'" _col(15) "bs`j'" _col(30) "KM weight_stab"
	if !fileexists("VBR/subgp_`k'/KM weight_stab/KM weight_stab subgp_`k' bs`j'.dta") {
	use "VBR/subgp_`k'/weight/weight subgp_`k' bs`j'.dta", clear
	count if invalid == 1
	if r(N) == 0 {
		stset tstop [pweight = weight_stab], enter(time tstart) failure(outcome)
		sts generate KM_s_w = s if arm_value == 1
		sts generate KM_ns_w = s if arm_value == 0
	}
	else {
		gen KM_s_w = .
		gen KM_ns_w = .
	}
	collapse (firstnm) KM_s_w KM_ns_w, by(tstop bs invalid)
	save "VBR/subgp_`k'/KM weight_stab/KM weight_stab subgp_`k' bs`j'.dta", replace
	}
}
}
*
* Finalize bootstrap datasets
qui forvalues k = 1/18 {
clear
forvalues j = 0/100 {
	append using "VBR/subgp_`k'/KM weight_stab/KM weight_stab subgp_`k' bs`j'.dta"
	*erase "unmatched/KM/KM subgp_`k' bs`j'.dta"
}
	gen hazard_s_w = 1 - KM_s_w
	gen hazard_ns_w = 1 - KM_ns_w
	gen odds_s_w = hazard_s_w/(1-hazard_s_w)
	gen odds_ns_w = hazard_ns_w/(1-hazard_ns_w)
	gen RR_w = hazard_s_w/hazard_ns_w
	gen diff_w = hazard_s_w - hazard_ns_w
	gen OR_w = odds_s_w / odds_ns_w
	rename tstop fup
	compress
	save "VBR/subgp_`k'/KM weight_stab subgp_`k' bs_all.dta", replace
}
*
cls
* KM estimate
qui forvalues k = 1/18 {
	use "VBR/subgp_`k'/KM weight_stab subgp_`k' bs_all.dta", clear
	
	keep if KM_s_w < . & KM_ns_w < .
	drop if invalid == 1
	bysort bs (fup) : keep if _n == _N
	keep if _n <= 1 + 100
	
	sum hazard_s_w if bs == 0, d
	scalar hazard_s_mean = r(mean)
	centile hazard_s_w if bs > 0, centile(2.5 97.5)
	scalar hazard_s_cil = r(c_1)
	scalar hazard_s_ciu = r(c_2)
	
	sum hazard_ns_w if bs == 0, d
	scalar hazard_ns_mean = r(mean)
	centile hazard_ns_w if bs > 0, centile(2.5 97.5)
	scalar hazard_ns_cil = r(c_1)
	scalar hazard_ns_ciu = r(c_2)
	
	noi di "subgp_`k'" _col(10) hazard_s_mean _col(25) hazard_s_cil _col(40) hazard_s_ciu _col(55) ///
	hazard_ns_mean _col(70) hazard_s_cil _col(85) hazard_s_ciu
}
*
* Absolute risk reduction
qui forvalues k = 1/18 {
	use "VBR/subgp_`k'/KM weight_stab subgp_`k' bs_all.dta", clear
	
	keep if KM_s_w < . & KM_ns_w < .
	drop if invalid == 1
	bysort bs (fup) : keep if _n == _N
	keep if _n <= 1 + 100
	
	sum diff_w if bs == 0, d
	scalar bs_mean = r(mean)
	scalar bs_p50 = r(p50)
	
	centile diff_w if bs > 0, centile(2.5 97.5)
	scalar bs_cil = r(c_1)
	scalar bs_ciu = r(c_2)
	noi di "subgp_`k'" _col(10) bs_mean _col(25) bs_p50 _col(40) bs_cil _col(55) bs_ciu
}
*
* N
qui forvalues k = 1/18 {
	use "VBR/subgp_`k'/cloned/cloned subgp_`k' bs0.dta", clear
	merge m:1 pseudo_key using "covid paxlovid VBR main.dta", keepusing(group) keep(1 3) nogen
	keep if (arm_value == 1 & group == 1) | (arm_value == 0 & group == 0)
	count if group == 1
	scalar N_1 = r(N)
	count if group == 0
	scalar N_0 = r(N)
	count if group == 1 & outcome == 1
	scalar n_e_1 = r(N)
	count if group == 0 & outcome == 1
	scalar n_e_0 = r(N)
	noi di "subgp_`k'" _col(15) N_1 _col(30) N_0 _col(45) n_e_1 _col(60) n_e_0
}
*
qui forvalues k = 1/18 {
	if `k' == 1 {
		noi di "subgp" _col(15) "N_1" _col(30) "N_0" _col(45) "n_e_1" _col(60) "n_e_0"
	}
	use "VBR/subgp_`k'/cloned/cloned subgp_`k' bs0.dta", clear
	merge m:1 pseudo_key using "covid paxlovid VBR IPTW subgp_`k'.dta", keep(3)
	count if arm_value == 1
	scalar N_1 = r(N)
	count if arm_value == 1 & outcome == 1
	scalar n_e_1 = r(N)
	count if arm_value == 0
	scalar N_0 = r(N)
	count if arm_value == 0 & outcome == 1
	scalar n_e_0 = r(N)
	
	noi di "subgp_`k'" _col(15) N_1 _col(30) N_0 _col(45) n_e_1 _col(60) n_e_0
}
*

cls
**# Viral burden rebound analysis - Sensitivity analyses
* sen1 - Day 0-2 vs Day >2
* VBR - sen1 - IPW
qui forvalues k = 1/18 {
	use "covid paxlovid VBR main.dta", clear
	merge m:1 pseudo_key using "covid paxlovid characteristics.dta", keep(3) nogen
	keep if subgp_`k' == 1

	replace group = inrange(time_to_treatment,0,2)
	logit group age gender charlson_index steroid_covid_bl immuno_hist covid_hist symptomatic_bl healthcare region_num_2 region_num_3 region_num_4 vaccine_status_2 vaccine_status_3 month_covid_cate_2 month_covid_cate_3 detection_cate inpatient_paxlovid censor_ctv_day_5

	predict prob_tx if e(sample)
	gen pscore_denom = prob_tx*group + (1-prob_tx)*(1-group)

	* generate numerator
	count
	scalar N_all = r(N)
	bysort group : gen N_group = _N
	gen pscore_num = N_group / N_all
	tab pscore_num group

	* generate stabilized IPTW
	gen _iptw_stab = pscore_num / pscore_denom
	sum _iptw_stab if group == 1, d
	noi di "subgp_`k'" _col(12) "gp=1" _col(20) "N="r(N) _col(30) r(mean) _col(45) r(sd) _col(60) r(min) _col(75) r(p25) _col(90) r(p50) _col(105) r(p90) _col(120) r(max)
	sum _iptw_stab if group == 0, d
	noi di "subgp_`k'" _col(12) "gp=0" _col(20) "N="r(N) _col(30) r(mean) _col(45) r(sd) _col(60) r(min) _col(75) r(p25) _col(90) r(p50) _col(105) r(p90) _col(120) r(max)
	sort pseudo_key

	compress
	save "covid paxlovid VBR IPTW sen1 subgp_`k'.dta", replace
}
*
* Target Trial Emulation
qui forvalues k = 1/18 {
	forvalues j = 0/200 {
		
	capture mkdir "VBR/sen1 subgp_`k'"
	capture mkdir "VBR/sen1 subgp_`k'/prepare"
	capture mkdir "VBR/sen1 subgp_`k'/cloned"
	capture mkdir "VBR/sen1 subgp_`k'/split"
	capture mkdir "VBR/sen1 subgp_`k'/Treatment"
	capture mkdir "VBR/sen1 subgp_`k'/NoTreatment"
	capture mkdir "VBR/sen1 subgp_`k'/weight"
	capture mkdir "VBR/sen1 subgp_`k'/KM weight_stab"
	
	*** Prepare dataset for analysis
	if `j'/10 == int(`j'/10) {
	noi di "subgp_`k'" _col(15) "bs`j'" _col(30) "setup"
	}
	* Setup for trial emulation
	if !fileexists("VBR/sen1 subgp_`k'/prepare/prepare subgp_`k' bs`j'.dta") {
	use "covid paxlovid VBR IPTW sen1 subgp_`k'.dta", clear
	
	* generate subgroups
	gen subgp_1 = 1

	gen subgp_2 = gender == 1
	gen subgp_3 = gender == 0

	gen subgp_4 = month_covid_cate == 1
	gen subgp_5 = month_covid_cate == 2
	gen subgp_6 = month_covid_cate == 3

	gen subgp_7 = vaccine_status > 1
	gen subgp_8 = vaccine_status == 1

	gen subgp_9 = inrange(charlson_index, 0, 6)
	gen subgp_10 = inrange(charlson_index, 7, 12)

	gen subgp_11 = date_steroid_covid <= date_index
	gen subgp_12 = !(date_steroid_covid <= date_index)

	gen subgp_13 = immuno_hist == 1
	gen subgp_14 = immuno_hist == 0

	gen subgp_15 = date_index == date_onset
	gen subgp_16 = date_index == date_covid
	
	gen subgp_17 = inpatient_paxlovid == 1
	gen subgp_18 = inpatient_paxlovid == 0

	keep if subgp_`k' == 1
	if `j' > 0 {
		set seed `j'
		bsample
	}
	gen date_event = date_VBR_end
	
	* VBR
	* organize / rename / generate variables
	gen time_to_event = date_event - date_index
	gen treatment = inrange(time_to_treatment, 0, 2)
	
	gen censor = 1 if time_to_treatment > 5
	gen date_censor = date_index + 5 if time_to_treatment > 5
	
	gen date_start_fu = date_index
	gen date_last_fu = min(date_death, mdy(02,12,2023), date_start_fu + 28, date_molnupiravir, date_censor)
	
	gen event = inrange(date_event, date_start_fu, date_last_fu)
	gen fup_obs = min(date_event-date_index, 28) if event == 1
	replace fup_obs = min(date_last_fu-date_index, 28) if event == 0

	* keep necessary variables
	keep pseudo_key group fup_obs event time_to_treatment time_to_event treatment censor
	gen bs = `j'
	compress
	save "VBR/sen1 subgp_`k'/prepare/prepare subgp_`k' bs`j'.dta", replace
	}
	
*** Cloning & censoring
	if `j'/10 == int(`j'/10) {
	noi di "subgp_`k'" _col(15) "bs`j'" _col(30) "cloning"
	}
	if !fileexists("VBR/sen1 subgp_`k'/cloned/cloned subgp_`k' bs`j'.dta") {
	* Prepare dataset for analysis
	use "VBR/sen1 subgp_`k'/prepare/prepare subgp_`k' bs`j'.dta", replace
	stset fup_obs, failure(event)

	* Arm A: no treatment within grace period (control: non-exposed group)
	gen outcomeA=_d // _d=event
	gen fupA=_t // _t=follow up time

	/// if the patient received treatment within day 0-1:
	/// 1. no event outcome, since the patient survived till censoring (treatment)
	replace outcomeA=0 if treatment==1 & time_to_treatment <= 2
	/// 2. follow up is censored at treatment
	replace fupA=time_to_treatment if treatment==1 & time_to_treatment <= 2

	* Arm B: treatment within grace period (treated: exposed group)
	gen outcomeB=_d 
	gen fupB=_t 

	/// if the patient survived the grace period and did not receive treatment within day 0-1:
	/// 1. no event outcome if the patient survived the grace period
	replace outcomeB=0 if (treatment==0 & _t>2) | (treatment==1 & time_to_treatment >2 & time_to_treatment !=.)
	/// 2. follow up is censored at end of grace period
	replace fupB=1 if (treatment==0 & _t>2) | (treatment==1 & time_to_treatment >2 & time_to_treatment != .)

	** append clones 
	preserve
		drop outcomeB fupB
		rename outcomeA outcome
		rename fupA fup
		gen arm="NoTreatment"
		tempfile a
		save "`a'", replace
	restore
		drop outcomeA fupA
		rename outcomeB outcome
		rename fupB fup
		gen arm="Treatment"	
		cap append using "`a'"

	// Weight models

	sort _all
	gen NewID=_n

	** add 1 day to 0-survivors
	*replace fup=1 if fup==0

	** exclude those initiated at day 0 in Late arm
	drop if time_to_treatment==0 & arm=="NoTreatment"

	** Weight model: define survival time and event indicator	
	* Early arm
	** Case 1: they do not deviate at time of treatment, but are not at risk of deviating any more
	gen wm_fup=time_to_treatment if arm=="Treatment" & time_to_treatment<=2 & time_to_treatment!=. & treatment==1 
	gen wm_outcome=0 if arm=="Treatment" & time_to_treatment<=2 & time_to_treatment!=. & treatment==1 

	** Case 2: they deviate at day 1
	replace wm_fup=2 if arm=="Treatment" & ((treatment==0 & fup>=2) | (time_to_treatment>2 & treatment==1))
	replace wm_outcome=1 if arm=="Treatment" & ((treatment==0 & fup>=2) | (time_to_treatment>2 & treatment==1)) 

	** Case 3: they do not deviate, but we need to keep their survival as observed and censor them as we do not know what happens appens afterwards
	replace wm_fup=fup if arm=="Treatment" & treatment==0 & fup<2
	replace wm_outcome=0 if arm=="Treatment" & treatment==0 & fup<2
	** add 1 days to 0-survivors
	replace wm_fup=1 if arm=="Treatment" & wm_fup==0 

	* Late arm
	** Case 1: they do not deviate at time of treatment, but are not at risk of deviating any more
	replace wm_fup=time_to_treatment if arm=="NoTreatment" & time_to_treatment<=2 & treatment==1 
	replace wm_outcome=1 if arm=="NoTreatment" & time_to_treatment<=2 & treatment==1 

	** Case 2: they deviate at day 1
	replace wm_fup=2 if arm=="NoTreatment" & ((treatment==0 & fup>=2) | (time_to_treatment>2 & treatment==1)) 
	replace wm_outcome=0 if arm=="NoTreatment" & ((treatment==0 & fup>=2) | (time_to_treatment>2 & treatment==1)) 

	** Case 3: they do not deviate, but we need to keep their survival as observed and censor them as we do not know what happens aerwards
	replace wm_fup=fup if arm=="NoTreatment" & treatment==0 & fup<2
	replace wm_outcome=0 if arm=="NoTreatment" & treatment==0 & fup<2

	** add 1 days to 0-survivors
	replace wm_fup=1 if arm=="NoTreatment" & wm_fup==0

	** censor patient at day 5 if initiated after day 5
	replace wm_outcome = 1 if time_to_treatment > 5
	
	gen arm_value = 1 if arm == "Treatment"
	replace arm_value = 0 if arm == "NoTreatment"
	drop arm
	order pseudo_key Early arm_value
	compress
	save "VBR/sen1 subgp_`k'/cloned/cloned subgp_`k' bs`j'.dta", replace
	}
	
* Split times
	if `j'/10 == int(`j'/10) {
	noi di "subgp_`k'" _col(15) "bs`j'" _col(30) "split"
	}
	if !fileexists("VBR/sen1 subgp_`k'/split/split subgp_`k' bs`j'.dta") {
	use "VBR/sen1 subgp_`k'/cloned/cloned subgp_`k' bs`j'.dta", clear
	
	** times of event
	stset fup, fail(outcome) id(NewID)
	stsplit, at(failures)

	gen tstart = _t0
	gen tstop = _t

	** times of censoring
	gen TrialEmul_cens = 1-outcome
	stset fup, fail(TrialEmul_cens) id(NewID)
	stsplit, at(failures)

	replace tstart = _t0 if tstart<_t0 & _t0 != . & _t != .
	replace tstop = _t if tstop>_t & _t0 != . & _t != .

	order pseudo_key arm tstart tstop
	sort _all
	compress
	save "VBR/sen1 subgp_`k'/split/split subgp_`k' bs`j'.dta", replace
	}
	
* IPCW - Early arm
	if `j'/10 == int(`j'/10) {
	noi di "subgp_`k'" _col(15) "bs`j'" _col(30) "Treatment"
	}
	if !fileexists("VBR/sen1 subgp_`k'/Treatment/Treatment subgp_`k' bs`j'.dta") {
	use "VBR/sen1 subgp_`k'/split/split subgp_`k' bs`j'.dta", clear
	keep if arm_value == 1
	merge m:1 pseudo_key using "covid paxlovid characteristics.dta", keep(3) nogen
	replace outcome = . if TrialEmul_cens == .

	* adapt records to the long format
	sort NewID tstart
	bysort NewID (tstart): replace wm_outcome = 0 if _n!=_N
	
	* Weight model:
	stset tstop, origin(time tstart) failure(wm_outcome) id(NewID)
	capture stcox age gender charlson_index steroid_covid_bl immuno_hist covid_hist symptomatic_bl healthcare region_num_2 region_num_3 region_num_4 vaccine_status_2 vaccine_status_3 month_covid_cate_2 month_covid_cate_3 detection_cate inpatient_paxlovid censor_ctv_day_5, efron
	if _rc == 0 {
		predict ch2_2, basech
		gen s0_2=exp(-ch2_2)
		predict xbCox2_2, xb
		gen weight = 1/s0_2^(exp(xbCox2_2)) 
		drop _st-xbCox2_2 
		gen invalid = 0
	}
	else {
		gen invalid = 1
		gen weight = .
		noi di "invalid"
	}
	keep pseudo_key arm tstart tstop event fup_obs time_to_treatment treatment bs outcome fup NewID wm_fup wm_outcome TrialEmul_cens weight invalid
	compress
	save "VBR/sen1 subgp_`k'/Treatment/Treatment subgp_`k' bs`j'.dta", replace
}

* IPCW - Late arm
	if `j'/10 == int(`j'/10) {
	noi di "subgp_`k'" _col(15) "bs`j'" _col(30) "NoTreatment"
	}
	if !fileexists("VBR/sen1 subgp_`k'/NoTreatment/Notreatment subgp_`k' bs`j'.dta") {
	use "VBR/sen1 subgp_`k'/split/split subgp_`k' bs`j'.dta", clear
	keep if arm_value == 0
	merge m:1 pseudo_key using "covid paxlovid characteristics.dta", keep(3) nogen
	replace outcome = . if TrialEmul_cens == .

	* adapt records to the long format
	sort NewID tstart
	bysort NewID: replace wm_outcome = 0 if _n!=_N

	* Weight model:
	stset tstop, origin(time tstart) failure(wm_outcome) id(NewID)
	capture stcox age gender charlson_index steroid_covid_bl immuno_hist covid_hist symptomatic_bl healthcare region_num_2 region_num_3 region_num_4 vaccine_status_2 vaccine_status_3 month_covid_cate_2 month_covid_cate_3 detection_cate inpatient_paxlovid censor_ctv_day_5, efron
	if _rc == 0 {
		predict ch2_2, basech
		gen s0_2=exp(-ch2_2)
		predict xbCox2_2, xb
		gen weight = 1/s0_2^(exp(xbCox2_2)) 
		drop _st-xbCox2_2 
		gen invalid = 0
	}
	else {
		gen invalid = 1
		gen weight = .
		noi di "invalid"
	}
	keep pseudo_key arm tstart tstop event fup_obs time_to_treatment treatment bs outcome fup NewID wm_fup wm_outcome TrialEmul_cens weight invalid
	compress
	save "VBR/sen1 subgp_`k'/NoTreatment/Notreatment subgp_`k' bs`j'.dta", replace
	}
	
* Combine & Generate weights
	if `j'/10 == int(`j'/10) {
	noi di "subgp_`k'" _col(15) "bs`j'" _col(30) "Combine"
	}
	if !fileexists("VBR/sen1 subgp_`k'/weight/weight subgp_`k' bs`j'.dta") {
	use "VBR/sen1 subgp_`k'/Treatment/Treatment subgp_`k' bs`j'.dta", clear
	append using "VBR/sen1 subgp_`k'/NoTreatment/Notreatment subgp_`k' bs`j'.dta"
	
	// create a new ID variable for each clone in each arm
	tostring NewID, gen(NewID_str)
	replace NewID_str = "00000" + NewID_str if length(NewID_str)== 1
	replace NewID_str = "0000" + NewID_str if length(NewID_str)== 2
	replace NewID_str = "000" + NewID_str if length(NewID_str)== 3
	replace NewID_str = "00" + NewID_str if length(NewID_str)== 4
	replace NewID_str = "0" + NewID_str if length(NewID_str)== 5
	gen Anal_ID = "1" + NewID_str if arm_value == 1
	replace Anal_ID = "2" + NewID_str if arm_value == 0

	replace weight = 1 if wm_outcome == 1 & tstop == 1
	sum weight
	local max_weight = r(max)

	keep pseudo_key arm_value tstart tstop treatment bs outcome wm_outcome Anal_ID weight invalid
	destring Anal_ID, replace

	* generate denominator
	rename weight _ipcw
	gen prob_denom = 1/_ipcw
	* generate numerator
	bysort arm_value tstart : gen N_group = _N
	bysort arm_value tstart : egen N_censored = sum(wm_outcome)
	gen prob_uncensored = (N_group-N_censored)/N_group
	gen prob_num = prob_uncensored
	sort arm_value pseudo_key tstart
	by arm_value pseudo_key : replace prob_num=prob_num*prob_num[_n-1] if _n!=1
	* generate stabilized IPCW
	gen _ipcw_stab = prob_num / prob_denom
	* generate combined weight
	merge m:1 pseudo_key using "covid paxlovid VBR IPTW sen1 subgp_`k'.dta", keep(3) keepusing(_iptw _iptw_stab)
	gen weight = _iptw * _ipcw
	gen weight_stab = _iptw_stab * _ipcw_stab
	compress
	save "VBR/sen1 subgp_`k'/weight/weight subgp_`k' bs`j'.dta", replace
	}
	
* Generate KM estimate
	if `j'/10 == int(`j'/10) {
	noi di "subgp_`k'" _col(15) "bs`j'" _col(30) "KM weight_stab"
	}
	if !fileexists("VBR/sen1 subgp_`k'/KM weight_stab/KM weight_stab subgp_`k' bs`j'.dta") {
	use "VBR/sen1 subgp_`k'/weight/weight subgp_`k' bs`j'.dta", clear
	count if invalid == 1
	if r(N) == 0 {
		stset tstop [pweight = weight_stab], enter(time tstart) failure(outcome)
		sts generate KM_s_w = s if arm_value == 1
		sts generate KM_ns_w = s if arm_value == 0
	}
	else {
		gen KM_s_w = .
		gen KM_ns_w = .
	}
	collapse (firstnm) KM_s_w KM_ns_w, by(tstop bs invalid)
	save "VBR/sen1 subgp_`k'/KM weight_stab/KM weight_stab subgp_`k' bs`j'.dta", replace
	}
}
}
*
* Finalize bootstrap datasets
qui forvalues k = 1/18 {
clear
forvalues j = 0/200 {
	capture append using "VBR/sen1 subgp_`k'/KM weight_stab/KM weight_stab subgp_`k' bs`j'.dta"
	*erase "unmatched/KM/KM subgp_`k' bs`j'.dta"
}
	gen hazard_s_w = 1 - KM_s_w
	gen hazard_ns_w = 1 - KM_ns_w
	gen odds_s_w = hazard_s_w/(1-hazard_s_w)
	gen odds_ns_w = hazard_ns_w/(1-hazard_ns_w)
	gen RR_w = hazard_s_w/hazard_ns_w
	gen diff_w = hazard_s_w - hazard_ns_w
	gen OR_w = odds_s_w / odds_ns_w
	rename tstop fup
	compress
	save "VBR/sen1 subgp_`k'/KM weight_stab subgp_`k' bs_all.dta", replace
}
*
* number of invalid
qui forvalues k = 1/18 {
	if `k' == 1 {
		noi di "subgp_`k'" _col(15) "n_all" _col(30) "n_valid" _col(45) "n_invalid" _col(60) "n_add"
	}
	use "VBR/sen1 subgp_`k'/KM weight_stab subgp_`k' bs_all.dta", clear
	keep if (KM_s_w < . & KM_ns_w < .) | invalid == 1
	bysort bs (fup) : keep if _n == _N
	count
	scalar n_all = r(N)
	count if invalid == 1
	scalar n_invalid = r(N)
	count if invalid == 0
	scalar n_valid = r(N)
	scalar n_add = max(101-n_valid, 0)
	noi di "subgp_`k'" _col(15) n_all _col(30) n_valid _col(45) n_invalid _col(60) n_add
}
*

cls
* KM estimate
qui forvalues k = 1/18 {
	use "VBR/sen1 subgp_`k'/KM weight_stab subgp_`k' bs_all.dta", clear
	
	keep if KM_s_w < . & KM_ns_w < .
	drop if invalid == 1
	bysort bs (fup) : keep if _n == _N
	keep if _n <= 1 + 100
	
	sum hazard_s_w if bs == 0, d
	scalar hazard_s_mean = r(mean)
	centile hazard_s_w if bs > 0, centile(2.5 97.5)
	scalar hazard_s_cil = r(c_1)
	scalar hazard_s_ciu = r(c_2)
	
	sum hazard_ns_w if bs == 0, d
	scalar hazard_ns_mean = r(mean)
	centile hazard_ns_w if bs > 0, centile(2.5 97.5)
	scalar hazard_ns_cil = r(c_1)
	scalar hazard_ns_ciu = r(c_2)
	
	noi di "subgp_`k'" _col(10) hazard_s_mean _col(25) hazard_s_cil _col(40) hazard_s_ciu _col(55) ///
	hazard_ns_mean _col(70) hazard_s_cil _col(85) hazard_s_ciu
}
*
* Absolute risk reduction
qui forvalues k = 1/18 {
	use "VBR/sen1 subgp_`k'/KM weight_stab subgp_`k' bs_all.dta", clear
	
	keep if KM_s_w < . & KM_ns_w < .
	drop if invalid == 1
	bysort bs (fup) : keep if _n == _N
	keep if _n <= 1 + 100
	
	sum diff_w if bs == 0, d
	scalar bs_mean = r(mean)
	scalar bs_p50 = r(p50)
	
	centile diff_w if bs > 0, centile(2.5 97.5)
	scalar bs_cil = r(c_1)
	scalar bs_ciu = r(c_2)
	noi di "subgp_`k'" _col(10) bs_mean _col(25) bs_p50 _col(40) bs_cil _col(55) bs_ciu
}
*
* N
qui forvalues k = 1/18 {
	use "VBR/sen1 subgp_`k'/cloned/cloned subgp_`k' bs0.dta", clear
	merge m:1 pseudo_key using "covid paxlovid VBR main.dta", keepusing(group) keep(1 3) nogen
	keep if (arm_value == 1 & group == 1) | (arm_value == 0 & group == 0)
	count if group == 1
	scalar N_1 = r(N)
	count if group == 0
	scalar N_0 = r(N)
	count if group == 1 & outcome == 1
	scalar n_e_1 = r(N)
	count if group == 0 & outcome == 1
	scalar n_e_0 = r(N)
	noi di "subgp_`k'" _col(15) N_1 _col(30) N_0 _col(45) n_e_1 _col(60) n_e_0
}
*
qui forvalues k = 1/18 {
	if `k' == 1 {
		noi di "subgp" _col(15) "N_1" _col(30) "N_0" _col(45) "n_e_1" _col(60) "n_e_0"
	}
	use "VBR/sen1 subgp_`k'/cloned/cloned subgp_`k' bs0.dta", clear
	merge m:1 pseudo_key using "covid paxlovid VBR IPTW sen1 subgp_`k'.dta", keep(3)
	count if arm_value == 1
	scalar N_1 = r(N)
	count if arm_value == 1 & outcome == 1
	scalar n_e_1 = r(N)
	count if arm_value == 0
	scalar N_0 = r(N)
	count if arm_value == 0 & outcome == 1
	scalar n_e_0 = r(N)
	
	noi di "subgp_`k'" _col(15) N_1 _col(30) N_0 _col(45) n_e_1 _col(60) n_e_0
}
*

cls
* sen2 - Day 0-3 vs Day >3
* VBR - sen2 - IPW
qui forvalues k = 1/18 {
use "covid paxlovid VBR main.dta", clear
merge m:1 pseudo_key using "covid paxlovid characteristics.dta", keep(3) nogen

keep if subgp_`k' == 1
replace group = inrange(time_to_treatment,0,3)
if `k' != 18 {
logit group age gender charlson_index steroid_covid_bl immuno_hist covid_hist symptomatic_bl healthcare region_num_2 region_num_3 region_num_4 vaccine_status_2 vaccine_status_3 month_covid_cate_2 month_covid_cate_3 detection_cate inpatient_paxlovid censor_ctv_day_5
}

if `k' == 18 {
	logit group gender charlson_index immuno_hist healthcare region_num_3 month_covid_cate_2
}
predict prob_tx if e(sample)
gen pscore_denom = prob_tx*group + (1-prob_tx)*(1-group)

* generate numerator
count
scalar N_all = r(N)
bysort group : gen N_group = _N
gen pscore_num = N_group / N_all
tab pscore_num group

* generate stabilized IPTW
gen _iptw_stab = pscore_num / pscore_denom
sum _iptw_stab if group == 1, d
noi di "subgp_`k'" _col(12) "gp=1" _col(20) "N="r(N) _col(30) r(mean) _col(45) r(sd) _col(60) r(min) _col(75) r(p25) _col(90) r(p50) _col(105) r(p90) _col(120) r(max)
sum _iptw_stab if group == 0, d
noi di "subgp_`k'" _col(12) "gp=0" _col(20) "N="r(N) _col(30) r(mean) _col(45) r(sd) _col(60) r(min) _col(75) r(p25) _col(90) r(p50) _col(105) r(p90) _col(120) r(max)
sort pseudo_key

compress
save "covid paxlovid VBR IPTW sen2 subgp_`k'.dta", replace
}
*
* Target Trial Emulation
qui forvalues k = 1/18 {
	forvalues j = 0/200 {
		
	capture mkdir "VBR/sen2 subgp_`k'"
	capture mkdir "VBR/sen2 subgp_`k'/prepare"
	capture mkdir "VBR/sen2 subgp_`k'/cloned"
	capture mkdir "VBR/sen2 subgp_`k'/split"
	capture mkdir "VBR/sen2 subgp_`k'/Treatment"
	capture mkdir "VBR/sen2 subgp_`k'/NoTreatment"
	capture mkdir "VBR/sen2 subgp_`k'/weight"
	capture mkdir "VBR/sen2 subgp_`k'/KM weight_stab"
	
	*** Prepare dataset for analysis
	if `j'/10 == int(`j'/10) {
	noi di "subgp_`k'" _col(15) "bs`j'" _col(30) "setup"
	}
	* Setup for trial emulation
	if !fileexists("VBR/sen2 subgp_`k'/prepare/prepare subgp_`k' bs`j'.dta") {
	use "covid paxlovid VBR IPTW sen2 subgp_`k'.dta", clear
	
	* generate subgroups
	gen subgp_1 = 1

	gen subgp_2 = gender == 1
	gen subgp_3 = gender == 0

	gen subgp_4 = month_covid_cate == 1
	gen subgp_5 = month_covid_cate == 2
	gen subgp_6 = month_covid_cate == 3

	gen subgp_7 = vaccine_status > 1
	gen subgp_8 = vaccine_status == 1

	gen subgp_9 = inrange(charlson_index, 0, 6)
	gen subgp_10 = inrange(charlson_index, 7, 12)

	gen subgp_11 = date_steroid_covid <= date_index
	gen subgp_12 = !(date_steroid_covid <= date_index)

	gen subgp_13 = immuno_hist == 1
	gen subgp_14 = immuno_hist == 0

	gen subgp_15 = date_index == date_onset
	gen subgp_16 = date_index == date_covid
	
	gen subgp_17 = inpatient_paxlovid == 1
	gen subgp_18 = inpatient_paxlovid == 0

	keep if subgp_`k' == 1
	sort pseudo_key
	if `j' > 0 {
		set seed `j'
		bsample
	}
	gen date_event = date_VBR_end
	
	* VBR
	* organize / rename / generate variables
	gen time_to_event = date_event - date_index
	gen treatment = inrange(time_to_treatment, 0, 2)
	
	gen censor = 1 if time_to_treatment > 5
	gen date_censor = date_index + 5 if time_to_treatment > 5
	
	gen date_start_fu = date_index
	gen date_last_fu = min(date_death, mdy(02,12,2023), date_start_fu + 28, date_molnupiravir, date_censor)
	
	gen event = inrange(date_event, date_start_fu, date_last_fu)
	gen fup_obs = min(date_event-date_index, 28) if event == 1
	replace fup_obs = min(date_last_fu-date_index, 28) if event == 0

	* keep necessary variables
	keep pseudo_key group fup_obs event time_to_treatment time_to_event treatment censor
	gen bs = `j'
	compress
	save "VBR/sen2 subgp_`k'/prepare/prepare subgp_`k' bs`j'.dta", replace
	}
	
	*** Cloning & censoring
	if `j'/10 == int(`j'/10) {
	noi di "subgp_`k'" _col(15) "bs`j'" _col(30) "cloning"
	}
	if !fileexists("VBR/sen2 subgp_`k'/cloned/cloned subgp_`k' bs`j'.dta") {
	* Prepare dataset for analysis
	use "VBR/sen2 subgp_`k'/prepare/prepare subgp_`k' bs`j'.dta", replace
	stset fup_obs, failure(event)

	* Arm A: no treatment within grace period (control: non-exposed group)
	gen outcomeA=_d // _d=event
	gen fupA=_t // _t=follow up time

	/// if the patient received treatment within day 0-3:
	/// 1. no event outcome, since the patient survived till censoring (treatment)
	replace outcomeA=0 if treatment==1 & time_to_treatment <= 3
	/// 2. follow up is censored at treatment
	replace fupA=time_to_treatment if treatment==1 & time_to_treatment <= 3

	* Arm B: treatment within grace period (treated: exposed group)
	gen outcomeB=_d 
	gen fupB=_t 

	/// if the patient survived the grace period and did not receive treatment within day 0-3:
	/// 1. no event outcome if the patient survived the grace period
	replace outcomeB=0 if (treatment==0 & _t>3) | (treatment==1 & time_to_treatment >3 & time_to_treatment !=.)
	/// 2. follow up is censored at end of grace period
	replace fupB=1 if (treatment==0 & _t>3) | (treatment==1 & time_to_treatment >3 & time_to_treatment != .)

	** append clones 
	preserve
		drop outcomeB fupB
		rename outcomeA outcome
		rename fupA fup
		gen arm="NoTreatment"
		tempfile a
		save "`a'", replace
	restore
		drop outcomeA fupA
		rename outcomeB outcome
		rename fupB fup
		gen arm="Treatment"	
		cap append using "`a'"

	// Weight models

	sort _all
	gen NewID=_n

	** exclude those initiated at day 0 in Late arm
	drop if time_to_treatment==0 & arm=="NoTreatment"

	** Weight model: define survival time and event indicator	
	* Early arm
	** Case 1: they do not deviate at time of treatment, but are not at risk of deviating any more
	gen wm_fup=time_to_treatment if arm=="Treatment" & time_to_treatment<=3 & time_to_treatment!=. & treatment==1 
	gen wm_outcome=0 if arm=="Treatment" & time_to_treatment<=3 & time_to_treatment!=. & treatment==1 

	** Case 2: they deviate at day 3
	replace wm_fup=3 if arm=="Treatment" & ((treatment==0 & fup>=3) | (time_to_treatment>3 & treatment==1))
	replace wm_outcome=1 if arm=="Treatment" & ((treatment==0 & fup>=3) | (time_to_treatment>3 & treatment==1)) 

	** Case 3: they do not deviate, but we need to keep their survival as observed and censor them as we do not know what happens appens afterwards
	replace wm_fup=fup if arm=="Treatment" & treatment==0 & fup<3
	replace wm_outcome=0 if arm=="Treatment" & treatment==0 & fup<3
	** add 1 days to 0-survivors
	replace wm_fup=1 if arm=="Treatment" & wm_fup==0 

	* Late arm
	** Case 1: they do not deviate at time of treatment, but are not at risk of deviating any more
	replace wm_fup=time_to_treatment if arm=="NoTreatment" & time_to_treatment<=3 & treatment==1 
	replace wm_outcome=1 if arm=="NoTreatment" & time_to_treatment<=3 & treatment==1 

	** Case 2: they deviate at day 3
	replace wm_fup=3 if arm=="NoTreatment" & ((treatment==0 & fup>=3) | (time_to_treatment>3 & treatment==1)) 
	replace wm_outcome=0 if arm=="NoTreatment" & ((treatment==0 & fup>=3) | (time_to_treatment>3 & treatment==1)) 

	** Case 3: they do not deviate, but we need to keep their survival as observed and censor them as we do not know what happens aerwards
	replace wm_fup=fup if arm=="NoTreatment" & treatment==0 & fup<3
	replace wm_outcome=0 if arm=="NoTreatment" & treatment==0 & fup<3

	** add 1 days to 0-survivors
	replace wm_fup=1 if arm=="NoTreatment" & wm_fup==0
	
	** censor patient at day 5 if initiated after day 5
	replace wm_outcome = 1 if time_to_treatment > 5
	
	gen arm_value = 1 if arm == "Treatment"
	replace arm_value = 0 if arm == "NoTreatment"
	drop arm
	order pseudo_key Early arm_value
	compress
	save "VBR/sen2 subgp_`k'/cloned/cloned subgp_`k' bs`j'.dta", replace
	}
	
	*** Split times
	if `j'/10 == int(`j'/10) {
	noi di "subgp_`k'" _col(15) "bs`j'" _col(30) "split"
	}
	if !fileexists("VBR/sen2 subgp_`k'/split/split subgp_`k' bs`j'.dta") {
	use "VBR/sen2 subgp_`k'/cloned/cloned subgp_`k' bs`j'.dta", clear
	
	** times of event
	stset fup, fail(outcome) id(NewID)
	stsplit, at(failures)

	gen tstart = _t0
	gen tstop = _t

	** times of censoring
	gen TrialEmul_cens = 1-outcome
	stset fup, fail(TrialEmul_cens) id(NewID)
	stsplit, at(failures)

	replace tstart = _t0 if tstart<_t0 & _t0 != . & _t != .
	replace tstop = _t if tstop>_t & _t0 != . & _t != .

	order pseudo_key arm tstart tstop
	sort _all
	compress
	save "VBR/sen2 subgp_`k'/split/split subgp_`k' bs`j'.dta", replace
	}
	
	*** IPCW - Early arm
	if `j'/10 == int(`j'/10) {
	noi di "subgp_`k'" _col(15) "bs`j'" _col(30) "Treatment"
	}
	if !fileexists("VBR/sen2 subgp_`k'/Treatment/Treatment subgp_`k' bs`j'.dta") {
	use "VBR/sen2 subgp_`k'/split/split subgp_`k' bs`j'.dta", clear
	keep if arm_value == 1
	merge m:1 pseudo_key using "covid paxlovid characteristics.dta", keep(3) nogen
	replace outcome = . if TrialEmul_cens == .

	* adapt records to the long format
	sort NewID tstart
	bysort NewID (tstart): replace wm_outcome = 0 if _n!=_N
	
	* Weight model:
	stset tstop, origin(time tstart) failure(wm_outcome) id(NewID)
	capture stcox age gender charlson_index steroid_covid_bl immuno_hist covid_hist symptomatic_bl healthcare region_num_2 region_num_3 region_num_4 vaccine_status_2 vaccine_status_3 month_covid_cate_2 month_covid_cate_3 detection_cate inpatient_paxlovid censor_ctv_day_5, efron
	if _rc == 0 {
		predict ch2_2, basech
		gen s0_2=exp(-ch2_2)
		predict xbCox2_2, xb
		gen weight = 1/s0_2^(exp(xbCox2_2)) 
		drop _st-xbCox2_2 
		gen invalid = 0
	}
	else {
		gen invalid = 1
		gen weight = .
		noi di "invalid"
	}
	keep pseudo_key arm tstart tstop event fup_obs time_to_treatment treatment bs outcome fup NewID wm_fup wm_outcome TrialEmul_cens weight invalid
	compress
	save "VBR/sen2 subgp_`k'/Treatment/Treatment subgp_`k' bs`j'.dta", replace
}

	*** IPCW - Late arm
	if `j'/10 == int(`j'/10) {
	noi di "subgp_`k'" _col(15) "bs`j'" _col(30) "NoTreatment"
	}
	if !fileexists("VBR/sen2 subgp_`k'/NoTreatment/Notreatment subgp_`k' bs`j'.dta") {
	use "VBR/sen2 subgp_`k'/split/split subgp_`k' bs`j'.dta", clear
	keep if arm_value == 0
	merge m:1 pseudo_key using "covid paxlovid characteristics.dta", keep(3) nogen
	replace outcome = . if TrialEmul_cens == .

	* adapt records to the long format
	sort NewID tstart
	bysort NewID: replace wm_outcome = 0 if _n!=_N

	* Weight model:
	stset tstop, origin(time tstart) failure(wm_outcome) id(NewID)
	capture stcox age gender charlson_index steroid_covid_bl immuno_hist covid_hist symptomatic_bl healthcare region_num_2 region_num_3 region_num_4 vaccine_status_2 vaccine_status_3 month_covid_cate_2 month_covid_cate_3 detection_cate inpatient_paxlovid censor_ctv_day_5, efron
	if _rc == 0 {
		predict ch2_2, basech
		gen s0_2=exp(-ch2_2)
		predict xbCox2_2, xb
		gen weight = 1/s0_2^(exp(xbCox2_2)) 
		drop _st-xbCox2_2 
		gen invalid = 0
	}
	else {
		gen invalid = 1
		gen weight = .
		noi di "invalid"
	}
	keep pseudo_key arm tstart tstop event fup_obs time_to_treatment treatment bs outcome fup NewID wm_fup wm_outcome TrialEmul_cens weight invalid
	compress
	save "VBR/sen2 subgp_`k'/NoTreatment/Notreatment subgp_`k' bs`j'.dta", replace
	}
	
	*** Combine & Generate weights
	if `j'/10 == int(`j'/10) {
	noi di "subgp_`k'" _col(15) "bs`j'" _col(30) "Combine"
	}
	if !fileexists("VBR/sen2 subgp_`k'/weight/weight subgp_`k' bs`j'.dta") {
	use "VBR/sen2 subgp_`k'/Treatment/Treatment subgp_`k' bs`j'.dta", clear
	append using "VBR/sen2 subgp_`k'/NoTreatment/Notreatment subgp_`k' bs`j'.dta"
	
	// create a new ID variable for each clone in each arm
	tostring NewID, gen(NewID_str)
	replace NewID_str = "00000" + NewID_str if length(NewID_str)== 1
	replace NewID_str = "0000" + NewID_str if length(NewID_str)== 2
	replace NewID_str = "000" + NewID_str if length(NewID_str)== 3
	replace NewID_str = "00" + NewID_str if length(NewID_str)== 4
	replace NewID_str = "0" + NewID_str if length(NewID_str)== 5
	gen Anal_ID = "1" + NewID_str if arm_value == 1
	replace Anal_ID = "2" + NewID_str if arm_value == 0

	replace weight = 1 if wm_outcome == 1 & tstop == 1
	sum weight
	local max_weight = r(max)

	keep pseudo_key arm_value tstart tstop treatment bs outcome wm_outcome Anal_ID weight invalid
	destring Anal_ID, replace

	* generate denominator
	rename weight _ipcw
	gen prob_denom = 1/_ipcw
	* generate numerator
	bysort arm_value tstart : gen N_group = _N
	bysort arm_value tstart : egen N_censored = sum(wm_outcome)
	gen prob_uncensored = (N_group-N_censored)/N_group
	gen prob_num = prob_uncensored
	sort arm_value pseudo_key tstart
	by arm_value pseudo_key : replace prob_num=prob_num*prob_num[_n-1] if _n!=1
	* generate stabilized IPCW
	gen _ipcw_stab = prob_num / prob_denom
	* generate combined weight
	merge m:1 pseudo_key using "covid paxlovid VBR IPTW sen2 subgp_`k'.dta", keep(3) keepusing(_iptw _iptw_stab)
	gen weight = _iptw * _ipcw
	gen weight_stab = _iptw_stab * _ipcw_stab
	compress
	save "VBR/sen2 subgp_`k'/weight/weight subgp_`k' bs`j'.dta", replace
	}
	
	*** Generate KM estimate
	if `j'/10 == int(`j'/10) {
	noi di "subgp_`k'" _col(15) "bs`j'" _col(30) "KM weight_stab"
	}
	if !fileexists("VBR/sen2 subgp_`k'/KM weight_stab/KM weight_stab subgp_`k' bs`j'.dta") {
	use "VBR/sen2 subgp_`k'/weight/weight subgp_`k' bs`j'.dta", clear
	count if invalid == 1
	if r(N) == 0 {
		stset tstop [pweight = weight_stab], enter(time tstart) failure(outcome)
		sts generate KM_s_w = s if arm_value == 1
		sts generate KM_ns_w = s if arm_value == 0
	}
	else {
		gen KM_s_w = .
		gen KM_ns_w = .
	}
	collapse (firstnm) KM_s_w KM_ns_w, by(tstop bs invalid)
	save "VBR/sen2 subgp_`k'/KM weight_stab/KM weight_stab subgp_`k' bs`j'.dta", replace
	}
}
}
*
* Finalize bootstrap datasets
qui forvalues k = 1/18 {
clear
forvalues j = 0/200 {
	capture append using "VBR/sen2 subgp_`k'/KM weight_stab/KM weight_stab subgp_`k' bs`j'.dta"
	*erase "KM/KM subgp_`k' bs`j'.dta"
}
	gen hazard_s_w = 1 - KM_s_w
	gen hazard_ns_w = 1 - KM_ns_w
	gen odds_s_w = hazard_s_w/(1-hazard_s_w)
	gen odds_ns_w = hazard_ns_w/(1-hazard_ns_w)
	gen RR_w = hazard_s_w/hazard_ns_w
	gen diff_w = hazard_s_w - hazard_ns_w
	gen OR_w = odds_s_w / odds_ns_w
	rename tstop fup
	compress
	save "VBR/sen2 subgp_`k'/KM weight_stab subgp_`k' bs_all.dta", replace
}
*
* number of invalid
qui forvalues k = 1/18 {
	if `k' == 1 {
		noi di "subgp_`k'" _col(15) "n_all" _col(30) "n_valid" _col(45) "n_invalid" _col(60) "n_add"
	}
	use "VBR/sen2 subgp_`k'/KM weight_stab subgp_`k' bs_all.dta", clear
	keep if (KM_s_w < . & KM_ns_w < .) | invalid == 1
	bysort bs (fup) : keep if _n == _N
	count
	scalar n_all = r(N)
	count if invalid == 1
	scalar n_invalid = r(N)
	count if invalid == 0
	scalar n_valid = r(N)
	scalar n_add = max(101-n_valid, 0)
	noi di "subgp_`k'" _col(15) n_all _col(30) n_valid _col(45) n_invalid _col(60) n_add
}
*
cls
* KM estimate
qui forvalues k = 1/18 {
	use "VBR/sen2 subgp_`k'/KM weight_stab subgp_`k' bs_all.dta", clear
	
	keep if KM_s_w < . & KM_ns_w < .
	drop if invalid == 1
	bysort bs (fup) : keep if _n == _N
	keep if _n <= 1 + 100
	
	sum hazard_s_w if bs == 0, d
	scalar hazard_s_mean = r(mean)
	centile hazard_s_w if bs > 0, centile(2.5 97.5)
	scalar hazard_s_cil = r(c_1)
	scalar hazard_s_ciu = r(c_2)
	
	sum hazard_ns_w if bs == 0, d
	scalar hazard_ns_mean = r(mean)
	centile hazard_ns_w if bs > 0, centile(2.5 97.5)
	scalar hazard_ns_cil = r(c_1)
	scalar hazard_ns_ciu = r(c_2)
	
	noi di "subgp_`k'" _col(10) hazard_s_mean _col(25) hazard_s_cil _col(40) hazard_s_ciu _col(55) ///
	hazard_ns_mean _col(70) hazard_s_cil _col(85) hazard_s_ciu
}
*
* Absolute risk reduction
qui forvalues k = 1/18 {
	use "VBR/sen2 subgp_`k'/KM weight_stab subgp_`k' bs_all.dta", clear
	
	keep if KM_s_w < . & KM_ns_w < .
	drop if invalid == 1
	bysort bs (fup) : keep if _n == _N
	keep if _n <= 1 + 100
	
	sum diff_w if bs == 0, d
	scalar bs_mean = r(mean)
	scalar bs_p50 = r(p50)
	
	centile diff_w if bs > 0, centile(2.5 97.5)
	scalar bs_cil = r(c_1)
	scalar bs_ciu = r(c_2)
	noi di "subgp_`k'" _col(10) bs_mean _col(25) bs_p50 _col(40) bs_cil _col(55) bs_ciu
}
*
* N
qui forvalues k = 1/18 {
	use "VBR/sen2 subgp_`k'/cloned/cloned subgp_`k' bs0.dta", clear
	merge m:1 pseudo_key using "covid paxlovid VBR main.dta", keepusing(group) keep(1 3) nogen
	keep if (arm_value == 1 & group == 1) | (arm_value == 0 & group == 0)
	count if group == 1
	scalar N_1 = r(N)
	count if group == 0
	scalar N_0 = r(N)
	count if group == 1 & outcome == 1
	scalar n_e_1 = r(N)
	count if group == 0 & outcome == 1
	scalar n_e_0 = r(N)
	noi di "subgp_`k'" _col(15) N_1 _col(30) N_0 _col(45) n_e_1 _col(60) n_e_0
}
*
qui forvalues k = 1/18 {
	if `k' == 1 {
		noi di "subgp" _col(15) "N_1" _col(30) "N_0" _col(45) "n_e_1" _col(60) "n_e_0"
	}
	use "VBR/sen2 subgp_`k'/cloned/cloned subgp_`k' bs0.dta", clear
	merge m:1 pseudo_key using "covid paxlovid VBR IPTW sen2 subgp_`k'.dta", keep(3)
	count if arm_value == 1
	scalar N_1 = r(N)
	count if arm_value == 1 & outcome == 1
	scalar n_e_1 = r(N)
	count if arm_value == 0
	scalar N_0 = r(N)
	count if arm_value == 0 & outcome == 1
	scalar n_e_0 = r(N)
	
	noi di "subgp_`k'" _col(15) N_1 _col(30) N_0 _col(45) n_e_1 _col(60) n_e_0
}
*
