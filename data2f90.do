/*stata file which reads cadastral maps from email of Aug 22 2020 and transform into input for fortran

clear all
cd "C:\Users\jbueren\Google Drive\overdrilling\fortran\Unobs_heterogeneity_ME\cadastral_maps"

*Select discretization of area types
local q=4

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
gen area_`q'_=1 if area_ac>=0.0 & area_ac<=1.3
replace area_`q'_=2 if area_ac>1.3 & area_ac<=2.3
replace area_`q'_=3 if area_ac>2.3 & area_ac<=4
replace area_`q'_=4 if area_ac>4 & area_ac<2000

keep pid map area_`q'_
reshape wide area_`q'_,i(pid) j(map)

foreach x of varlist area_`q'_1-area_`q'_14 {
  replace `x' = -9 if (`x' >= .)
}
sort pid
export delimited area_* using "fortran_files\area_type.csv",replace novarnames

*number of plots in each map
forval m=1/14{
count if area_4_`m'!=-9
}
*/
*Primitives on flow and failure by type and monsoon: I generate 4 different types of unobserved heterogeneity.
clear all
cd "C:\Users\jbueren\Google Drive\overdrilling\fortran\Unobs_heterogeneity_ME"
import excel using "primitives\flow_fail_prob_r.xls",firstrow   

egen id=group(N M T)

reshape long Pfail_,i(id) j(T2)
drop id
egen T3=group(T T2)
br
drop T T2
rename T3 T
order T,after(M)
br
export delimited using "primitives\flow_fail_prob_r",novarnames  replace 

*Primitives: probability of success
clear all
cd "C:\Users\jbueren\Google Drive\overdrilling\fortran\Unobs_heterogeneity_ME"
import excel using "data\drill_export_r.xls",firstrow   
collapse P_R P_S,by(map_village)
encode map_village,g(nb)
drop if nb==.
export delimited nb P_R P_S using "primitives\rain_success_pr.csv",replace novarnames nolabel 
br


*Estimation data
clear all
cd "C:\Users\jbueren\Google Drive\overdrilling\fortran\Unobs_heterogeneity_ME\data"
import excel using "drill_export_r.xls",firstrow   
encode map_village,g(nb)
drop if nb==.
br
*Number of well in the plot
gen n=ref1well+2*ref2well
replace drill=-9 if n==2

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
bys a_type: sum area,d

gen P_type=min(Nplots_adj,6)

sort RespondentID year

gen P_T1=Pflow_T1_e*Pfail1_T1_e
gen P_T2=Pflow_T1_e*Pfail2_T1_e
gen P_T3=Pflow_T2_e*Pfail1_T2_e
gen P_T4=Pflow_T2_e*Pfail2_T2_e

*replace drill=-9 if IN_FF_SMPL!=1

export delimited nb P_type a_type n f0_N - f10_N P_T1 P_T2 P_T3 P_T4 drill using "drill_export_r.csv",replace novarnames nolabel 

*statistics by area en number of wells around
bys a_type: sum drill if drill>=0

gen max_pr=max(f0_N, f1_N, f2_N, f3_N, f4_N, f5_N ,f6_N ,f7_N ,f8_N, f9_N ,f10_N)
gen big_N=.
forval i=0/10{
replace big_N=`i' if f`i'_N==max_pr
}
drop max_pr
bys big_N: sum drill if drill>=0

gen max_pr=max(P_T1, P_T2, P_T3, P_T4)
gen modal_type=.
forval i=1/4{
replace modal_type=`i' if P_T`i'==max_pr
}
drop max_pr
bys modal_type: sum drill if drill>=0 

bys a_type n: sum drill if drill>=0



