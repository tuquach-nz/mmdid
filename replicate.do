*** These codes are to replicate the results in the paper ""
*** Run mmdid.ado package
*** mmdid means Multilevel Modeling Difference-in-Differences estimator

run "mmdid.ado"

*** Install CSDID package if needed

ssc install csdid
ssc install drdid

/*
Replicate the results in section 6.1
After drawing the graph using the coefplot command, we continue to edit
the graphs using the graph editor in Stata.
The information in Table 2 is also harvested using the codes below.
*/
use "mnw.dta"

/* Define event time variable as follows

gen ei = cond(treat==1, year - first_treat, -8)

*/

*** Generate lag and lead variables in TWFE specification
*** Lag variables
forv i =0/3 {
gen pp`i' = cond(ei==`i',1,0)
}

*** Lead variables
forv i = 0/6 {
gen pm`i' = cond(ei==-`i',1,0)
}

*** Run TWFE specification and graph Fig 1.a
xtset countyreal year

xtreg lnemp18 treat lnpop pm6 pm5 pm4 pm3 pm2  pp0 pp1 pp2 pp3 i.year, fe cluster( countyreal)

nlcom (_b[pm6]) (_b[pm5]) (_b[pm4]) (_b[pm3]) (_b[pm2]) (0) (_b[pp0]) (_b[pp1]) (_b[pp2]) (_b[pp3]), post level(95)
estimate store twfe

coefplot (twfe, ciopt(recast(rcap)) recast(connected) level(95)), vertical keep(_nl_*) ///
			coeflabel(_nl_1=-6 _nl_2=-5 _nl_3=-4 _nl_4=-3 _nl_5=-2 ///
			_nl_6=-1 _nl_7=0 _nl_8=1 _nl_9=2 _nl_10=3)  xlabel(,labsize(med)) ///
			yline(0, lcolor(black) lwidth(thin) lpattern(dash)) ylabel(,labsize(med) ) ///
			graphregion(color(white)) legend(off) xtitle("Event time") ytitle("ATT")
			
graph export "twfe_emp1.pdf", as(pdf) name("Graph")

*** Run Callaway and Sant'Anna approach and graph Fig 1.b

csdid lnemp18 lnpop, i( countyreal) t( year) gvar( first_treat)
estat event, post
estimate store csdid

coefplot (csdid, ciopt(recast(rcap)) recast(connected) level(95)), vertical keep(Tm* Tp*) ///
			coeflabel(Tm5=-5 Tm4=-4 Tm3=-3 Tm2=-2 Tm1=-1 Tp0=0 Tp1=1 Tp2=2 Tp3=3) ///
			xlabel(,labsize(med)) yline(0, lcolor(black) lwidth(thin) lpattern(dash)) ///
			ylabel(,labsize(med) ) graphregion(color(white)) legend(off) xtitle("Event time") ytitle("ATT")
			
graph export "csdid_emp1.pdf", as(pdf) name("Graph")

*** Run mmdid approach and graph Fig 1.c (homogeneity) and 1.d (heterogeneity)
*** Generate the interaction variable between the independent variable and the treatment variable

gen cov = treat*lnpop
xi: mmdid lnemp18 lnpop cov i.year treat, tv( treat) ev(ei) 
estimate store mmdid

coefplot (mmdid, ciopt(recast(rcap)) recast(connected) level(95)), vertical ///
			xlabel(,labsize(med)) yline(0, lcolor(black) lwidth(thin) lpattern(dash)) ///
			ylabel(,labsize(med) ) graphregion(color(white)) legend(off) xtitle("Event time") ytitle("ATT")

graph export "mmdid_emp1_homo.pdf", as(pdf) name("Graph")

xi: mmdid lnemp18 lnpop cov i.year treat, tv( treat) ev(ei) att( first_treat)
estimate store mmdid

coefplot (mmdid, ciopt(recast(rcap)) recast(connected) level(95)), vertical ///
			xlabel(,labsize(med)) yline(0, lcolor(black) lwidth(thin) lpattern(dash)) ///
			ylabel(,labsize(med) ) graphregion(color(white)) legend(off) xtitle("Event time") ytitle("ATT")

graph export "mmdid_emp1_hete.pdf", as(pdf) name("Graph")

/*
Replicate the results in section 6.2

*/

use "hospitalization.dta"

*** Table 3

sum age_hosp male white black hispanic insured_pv if wave==7
tabstat oop_spend riearnsemp if wave <11, by(wave) stat( n mean sd)

*** We run the TWFE specification and the mmdid approach for the case where there is no control group.
*** Generate lag and lead variables in TWFE specification

*** Lag variables
forv i =0/3 {
gen pp`i' = cond(ei==`i',1,0)
}

*** Lead variables
forv i = 0/4 {
gen pm`i' = cond(ei==-`i',1,0)
}

*** Table 4 - the first two columns in panel B
*** TWFE specification
xtset hhidpn wave
xtreg  oop_spend pp3 pp2 pp1 pp0 pm2 pm3 i.wave, fe

*** mmdid approach
xi: mmdid oop_spend ever_hospitalized i.wave, tv(ever_hospitalized) ev(ei) att(start)

*** We used the individuals who started hospitalization in wave 11 as the control group. 
*** Therefore, we drop observations in wave 11 and modify the treatment variable.

preserve

replace ever_hospitalized =0 if start==11
drop if wave==11
replace ei =-5 if ever_hospitalized==0

*** The last two columns in panel B
xi: mmdid oop_spend ever_hospitalized i.wave, tv(ever_hospitalized) ev(ei) att(start)
csdid  oop_spend, i( hhidpn) t(wave) gvar(start)

*** Panel C in Table 4

xi: mmdid riearnsemp ever_hospitalized i.wave, tv(ever_hospitalized) ev(ei) att(start)
xi: mmdid riearnsemp ever_hospitalized i.wave, tv(ever_hospitalized) ev(ei) att(start insured_pv)

csdid riearnsemp if insured_pv == 0, i( hhidpn) t(wave) gvar(start)
csdid riearnsemp if insured_pv == 100, i( hhidpn) t(wave) gvar(start)

*** Figure 2

xi: mmdid oop_spend ever_hospitalized i.wave, tv(ever_hospitalized) ev(ei) att(start)
mmdid_plot, sep
graph name a

xi: mmdid riearnsemp ever_hospitalized i.wave, tv(ever_hospitalized) ev(ei) att(start)
mmdid_plot, sep
graph name b

graph combine a b

restore

*** End ***