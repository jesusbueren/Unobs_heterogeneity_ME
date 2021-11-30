*stata file which reads cadastral maps from email of Aug 22 2020 and transform into input for fortran

clear all
cd "C:\Users\jbueren\Google Drive\overdrilling\fortran\Unobs_heterogeneity_ME\cadastral_maps"

*Map
forval m=1/14{
clear all
import delimited using "nbs\nbs_`m'.csv"
drop if _n==1
drop v1
export delimited using "fortran_files\map_`m'.csv",replace novarnames
}

*Area types
forval m=1/14{
clear all
import delimited using "tbl\\`m'_attrib.csv"
sort pid
keep pid area_ac
gen map=`m'
if `m'>1{
append using `cal'
}
tempfile cal
save `cal'
}
drop if pid==.

rename area_ac area_ac_
keep pid map area_ac_
reshape wide area_ac_,i(pid) j(map)

foreach x of varlist area_ac_1-area_ac_14 {
  replace `x' = -9 if (`x' >= .)
}
sort pid
export delimited area_* using "fortran_files\area_type.csv",replace novarnames

*number of plots in each map
forval m=1/14{
count if area_ac_`m'!=-9
}
*/
*Primitives on flow and failure by type and monsoon: I generate 4 different types of unobserved heterogeneity.
clear all
cd "C:\Users\jbueren\Google Drive\overdrilling\fortran\Unobs_heterogeneity_ME"
*import excel using "primitives\allprobs_lnN_rho_median.xls",firstrow  
import excel using "primitives\allprobs_lnN_rho.xls",firstrow  
/*import excel using "primitives\flow_lnN_fail_N_prob_2T.xls",firstrow  
sort T N M
forval i=3/4{
set obs `=_N+18*2'
replace N=N[_n-18*2*2] if N==.
replace M=M[_n-18*2*2] if M==.
replace T=`i' if T==.
if `i'==3{
foreach var of varlist Pq1-Pq5 {
replace `var'=`var'[_n-18*2*2] if `var'==.
}
replace Pfail=Pfail[_n-18*2] if Pfail==.
}
if `i'==4{
foreach var of varlist Pq1-Pq5 {
replace `var'=`var'[_n-18*2*2] if `var'==.
}
replace Pfail=Pfail[_n-18*2*3] if Pfail==.
}
}
*/
sort N M T 
br
/* for old files
egen id=group(N M T)

reshape long Pfail_,i(id) j(T2)
drop id
egen T3=group(T T2)
br
drop T T2
rename T3 T
order T,after(M)
br
*/
export delimited using "primitives\flow_fail_prob_r",novarnames  replace 

*Primitives: probability of success
clear all
cd "C:\Users\jbueren\Google Drive\overdrilling\fortran\Unobs_heterogeneity_ME"
import excel using "data\drill_export_lnN_3T.xls",firstrow   
collapse P_R P_S,by(map_village)
encode map_village,g(nb)
drop if nb==.
export delimited nb P_R P_S using "primitives\rain_success_pr.csv",replace novarnames nolabel 
br


*Estimation data
clear all
cd "C:\Users\jbueren\Google Drive\overdrilling\fortran\Unobs_heterogeneity_ME\data"
import excel using "drill_export_lnN_3T.xls",firstrow   
encode map_village,g(nb)
drop if nb==.
br
*Number of well in the plot
gen n=ref1well+2*ref2well
replace drill=0 if n==2

sort RespondentID year
by RespondentID: gen delta_n=n-n[_n-1]
replace n=n-1 if delta_n==2
drop delta_n
by RespondentID: gen delta_n=n-n[_n-1]
tab delta_n if drill[_n-1]==0
by RespondentID:  replace drill=1 if delta_n[_n+1]==1 & drill==0

*Measurement error taking into account knowledge of number of well in own plot
gen f9=0
gen f10=0
forval i=0/10{
gen f`i'_N=0
}
forval i=0/10{
replace f`i'_N=f`i' if n==0
}
forval i=1/10{
local j=`i'-1
replace f`i'_N=f`j' if n==1
}
forval i=2/10{
local j=`i'-2
replace f`i'_N=f`j' if n==2
}

local q=4
xtile a_type=area,n(`q')


recode a_type (1/2=1)(3/4=2)

bys a_type: sum area,d

gen P_type=min(Nplots_adj,6)

sort RespondentID year

rename Pflow_T1 P_T1
rename Pflow_T2 P_T2
rename Pflow_T3 P_T3

by RespondentID: egen Total_n=total(n)
*drop if Total_n==10
/*
replace P_T1=1/3 if IMPUTE==1
replace P_T2=1/3 if IMPUTE==1
replace P_T3=1/3 if IMPUTE==1
*/

by RespondentID: egen total_attempts=total(drill) if drill!=-9
by RespondentID: egen total_attempts_2=mean(total_attempts)

replace IMPUTE=0
*replace IMPUTE=1 if total_attempts_2==0 & n==0

sort RespondentID year
gen can_be_zombie=0
by RespondentID: replace can_be_zombie=1 if total_attempts_2==0 & n[1]==0

*replace drill=0 if n==1 & a_type==1
*replace n=1 if n==2 & a_type==1
replace P_type=3  if P_type==2
export delimited nb P_type a_type n f0_N - f10_N P_T1 P_T2 P_T3 drill IMPUTE can_be_zombie using "drill_export_r.csv",replace novarnames nolabel 

*statistics by area en number of wells around
bys a_type: sum drill if drill>=0

gen max_pr=max(f0, f1, f2, f3, f4, f5 ,f6 ,f7 ,f8, f9 ,f10)
gen big_N=.
forval i=0/10{
replace big_N=`i' if f`i'_N==max_pr
}
drop max_pr
bys big_N: sum drill if drill>=0

bys a_type n: sum drill if drill>=0

bys RespondentID: egen N_bar=mean(big_N)
reg drill i.n big_N N_bar if can_be_zombie==0 & n==0

bys RespondentID: egen total_d=total(drill)
reg total_d N_bar n if year==2012 & can_be_zombie==0 

clear all
cd "C:\Users\jbueren\Google Drive\overdrilling\fortran\Unobs_heterogeneity_ME\Unobs_heterogeneity"
infix i_l 1-5 t_l 6-10 P 11-15 A 16-20 n_l 21-25 N 26-30 drill 31-35 using panel.txt
replace N=N-1
replace n_l=n_l-1
replace drill=. if drill==-9
sort i_l t_l
by i_l: gen delta_n=n[_n+1]-n
by i_l: egen total_d=total(drill)
by i_l: egen N_bar=mean(N)
reg drill i.n_l N N_bar 
reg total_d N_bar n if t_l==1

