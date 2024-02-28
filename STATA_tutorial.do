/*PROJECT: WHO PROJECT
AUTHOR: JOSEPH CHALLENGER
DATE: AUGUST 2023
*/
* Change working directory if needed
*e.g. 
*cd C:\Users\username\Documents

import delimited "example_dataset.csv"

tab tot_dead
summarize tot_dead
levelsof(hut)
levelsof(sleeper)
levelsof(treatment)

/* Determine the non-inferiority margin (NIM)
NIM: mosquito mortality induced by the candidate net should be no more than 7% less
 than that induced by the active comparator net.
First, calculate the (unadjusted) mosquito mortality in each trial arm */
collapse (sum) sum1=tot_dead sum2=total, by(treatment)
gen prop_dead = sum1/sum2
list
gen or1 = (prop_dead - 0.07)/(1-prop_dead + 0.07)
gen or2 = (prop_dead)/(1-prop_dead)
*Calculate the odds-ratio (OR) for the NIM
gen nim = or1/or2
list
save "aggregated_mortality.dta"

/* We will need alternative NIMs for the combined analyses 
(i.e. considering washed and unwashed nets together */
clear
import delimited "example_dataset.csv"
collapse (sum) sum1=tot_dead sum2=total, by(itn)
gen prop_dead = sum1/sum2
list
gen or1 = (prop_dead - 0.07)/(1-prop_dead + 0.07)
gen or2 = (prop_dead)/(1-prop_dead)
*Calculate the odds-ratio (OR) for the NIM
gen nim = or1/or2
list
save "aggregated_mortality_itn.dta"

/* We will also need NIM for the blood feeding analysis. 
In this case, blood feeding in the candidate arm should be no more than 7% 
more than in the active comparator arm*/
clear
import delimited "example_dataset.csv"
collapse (sum) sum1=tot_bf sum2=total, by(treatment)
gen prop_fed = sum1/sum2
list
gen or1 = (prop_fed + 0.07)/(1-prop_fed - 0.07)
gen or2 = (prop_fed)/(1-prop_fed)
*Calculate the odds-ratio (OR) for the NIM
gen nim = or1/or2
list
save "aggregated_bf.dta"

/* Finally, we will need NIM for the combined blood feeding analysis. */
clear
import delimited "example_dataset.csv"
collapse (sum) sum1=tot_bf sum2=total, by(itn)
gen prop_fed = sum1/sum2
list
gen or1 = (prop_fed + 0.07)/(1-prop_fed - 0.07)
gen or2 = (prop_fed)/(1-prop_fed)
*Calculate the odds-ratio (OR) for the NIM
gen nim = or1/or2
list
save "aggregated_bf_itn.dta"

/* ############################################################## 
			1. Mosquito mortality (unwashed ITNS)     
 ############################################################## */

/* A. Reload dataset to carry out the regression*/
clear
import delimited "example_dataset.csv"
append using "aggregated_mortality.dta"
*Remove variables we don't need anymore
drop sum1 sum2 prop_dead or1 or2

/* B. At the moment, 'treatment' is a string variable. We need a 
factor variable for the regression model */
encode(treatment), generate(treatment2)
replace treatment2=. if day==.
levelsof(treatment2)
*See how the treatment2 levels correspond to treatment*
label list treatment2

/* C. For the regression model we use the function 'blogit', 
which fits a logistic regression model to aggregated data */
blogit tot_dead total i.treatment2 i.hut i.sleeper i.day

/* In this case, we wish to use treatment2=1 as the baseline category, so this is OK.
 If we wished to change the baseline category for treatment (.e.g. to treatment2=6) 
 we could do this: blogit tot_dead total i.hut ib6.treatment2 i.sleeper i.day */

*Here is how the model is stored in Stata's memory*
ereturn list

*Calculate the odds ratio (OR) and 95% CI for the unwashed candidate net
gen or_model = exp(_b[_outcome:3.treatment2])
gen or_model_lower = exp(_b[_outcome:3.treatment2] - 1.96* _se[_outcome:3.treatment2])
gen or_model_upper = exp(_b[_outcome:3.treatment2] + 1.96* _se[_outcome:3.treatment2])

*Alternatively, we could have asked Stata to calculate the ORs for us, like this:*
blogit tot_dead total i.treatment2 i.hut i.sleeper i.day, or

*Recall the non-inferiority margin. Use this to compare with the OR calculated above
list if treatment=="Active_comparator_unwashed" & missing(day)

/* D. Now check whether the unwashed candidate net is superior to the unwashed
standard comparator.  */
blogit tot_dead total ib6.treatment2 i.hut i.sleeper i.day

/* ############################################################## 
			2. Mosquito mortality (washed ITNS)     
 ############################################################## */
 
/* A. Clear memory and reload data */
clear
import delimited "example_dataset.csv"
append using "aggregated_mortality.dta"
*Remove variables we don't need anymore
drop sum1 sum2 prop_dead or1 or2

/* B. At the moment, 'treatment' is a string variable. We need a 
factor variable for the regression model */
encode(treatment), generate(treatment2)
replace treatment2=. if day==.
levelsof(treatment2)
*See how the treatment2 levels correspond to treatment*
label list treatment2

/* C. For the regression model we use the function 'blogit', 
which fits a logistic regression model to aggregated data. 
We use level 2 of 'treatment2' as our baseline category
(this is 'Active_comparator_washed' */
blogit tot_dead total ib2.treatment2 i.hut i.sleeper i.day

/* Calculate the odds ratio (OR) and 95% CI for the unwashed candidate net,
 which is level 4 of treatment2 */
gen or_model = exp(_b[_outcome:4.treatment2])
gen or_model_lower = exp(_b[_outcome:4.treatment2] - 1.96* _se[_outcome:4.treatment2])
gen or_model_upper = exp(_b[_outcome:4.treatment2] + 1.96* _se[_outcome:4.treatment2])

*Alternatively, we could have asked Stata to calculate the ORs for us, like this:*
blogit tot_dead total ib2.treatment2 i.hut i.sleeper i.day, or

*Recall the non-inferiority margin. Use this to compare with the OR calculated above
list if treatment=="Active_comparator_washed" & missing(day)

/* D. Now check whether the washed candidate net is superior to the washed
standard comparator (level 7 of 'treatment2').  */
blogit tot_dead total ib7.treatment2 i.hut i.sleeper i.day

 
/* ############################################################## 
		3. Mosquito mortality (combined unwashed & washed ITNS)     
 ############################################################## */
 
/* A. Clear memory and reload data */
clear
import delimited "example_dataset.csv"
append using "aggregated_mortality_itn.dta"
*Remove variables we don't need anymore
drop sum1 sum2 prop_dead or1 or2
 
/* B. At the moment, 'itn' is a string variable. We need a 
factor variable for the regression model */
encode(itn), generate(itn2)
replace itn2=. if day==.
levelsof(itn2)
*See how the itn2 levels correspond to itn*
label list itn2

/* C. For the regression model we use the function 'blogit', 
which fits a logistic regression model to aggregated data. 
Use level 1 of itn2 for the baseline category (i.e. the default) */
blogit tot_dead total i.itn2 i.hut i.sleeper i.day i.wash

/* Calculate the odds ratio (OR) and 95% CI for the candidate net,
 which is level 2 of itn2 */
gen or_model = exp(_b[_outcome:2.itn2])
gen or_model_lower = exp(_b[_outcome:2.itn2] - 1.96* _se[_outcome:2.itn2])
gen or_model_upper = exp(_b[_outcome:2.itn2] + 1.96* _se[_outcome:2.itn2])

*Alternatively, we could have asked Stata to calculate the ORs for us, like this:*
blogit tot_dead total i.itn2 i.hut i.sleeper i.day i.wash, or

*Recall the non-inferiority margin. Use this to compare with the OR calculated above
list if itn=="Active_comparator" & missing(day)

/* D. Now check whether the candidate net is superior to the
standard comparator (level 4 of 'itn2').  */
blogit tot_dead total ib4.itn2 i.hut i.sleeper i.day i.wash
 
/* ############################################################## 
			4. Blood feeding (unwashed ITNS)     
 ############################################################## */
 
/* A. Clear memory and reload data */
clear
import delimited "example_dataset.csv"
append using "aggregated_bf.dta"
*Remove variables we don't need anymore
drop sum1 sum2 prop_fed or1 or2

/* B. At the moment, 'treatment' is a string variable. We need a 
factor variable for the regression model */
encode(treatment), generate(treatment2)
replace treatment2=. if day==.
levelsof(treatment2)
*See how the treatment2 levels correspond to treatment*
label list treatment2

/* C. For the regression model we use the function 'blogit', 
which fits a logistic regression model to aggregated data */
blogit tot_bf total i.treatment2 i.hut i.sleeper i.day

/* Calculate the odds ratio (OR) and 95% CI for the candidate net,
 which is level 3 of treatment2 */
gen or_model = exp(_b[_outcome:3.treatment2])
gen or_model_lower = exp(_b[_outcome:3.treatment2] - 1.96* _se[_outcome:3.treatment2])
gen or_model_upper = exp(_b[_outcome:3.treatment2] + 1.96* _se[_outcome:3.treatment2])

*Alternatively, we could have asked Stata to calculate the ORs for us, like this:*
blogit tot_bf total i.treatment2 i.hut i.sleeper i.day, or

*Recall the non-inferiority margin. Use this to compare with the OR calculated above
list if treatment=="Active_comparator_unwashed" & missing(day)

/* D. Now check whether the unwashed candidate net is superior to the unwashed
standard comparator (level 4 of 'itn2').  */
blogit tot_bf total ib6.treatment2 i.hut i.sleeper i.day

/* ############################################################## 
			5. Blood feeding  (washed ITNS)     
 ############################################################## */
 
/* A. Clear memory and reload data */
clear
import delimited "example_dataset.csv"
append using "aggregated_bf.dta"
*Remove variables we don't need anymore
drop sum1 sum2 prop_fed or1 or2

/* B. At the moment, 'treatment' is a string variable. We need a 
factor variable for the regression model */
encode(treatment), generate(treatment2)
replace treatment2=. if day==.
levelsof(treatment2)
*See how the treatment2 levels correspond to treatment*
label list treatment2

/* C. For the regression model we use the function 'blogit', 
which fits a logistic regression model to aggregated data. Set the
baseline category to 'Active_comparator_washed' */
blogit tot_bf total ib2.treatment2 i.hut i.sleeper i.day

/* Calculate the odds ratio (OR) and 95% CI for the candidate net,
 which is level 3 of treatment2 */
gen or_model = exp(_b[_outcome:4.treatment2])
gen or_model_lower = exp(_b[_outcome:4.treatment2] - 1.96* _se[_outcome:4.treatment2])
gen or_model_upper = exp(_b[_outcome:4.treatment2] + 1.96* _se[_outcome:4.treatment2])

*Alternatively, we could have asked Stata to calculate the ORs for us, like this:*
blogit tot_bf total ib2.treatment2 i.hut i.sleeper i.day, or

*Recall the non-inferiority margin. Use this to compare with the OR calculated above
list if treatment=="Active_comparator_washed" & missing(day)

/* D. Now check whether the washed candidate net is superior to the washed
standard comparator (level 4 of 'itn2').  */
blogit tot_bf total ib7.treatment2 i.hut i.sleeper i.day
 
/* ############################################################## 
		6. Blood feeding (combined unwashed & washed ITNS)     
############################################################## */

/* A. Clear memory and reload data */
clear
import delimited "example_dataset.csv"
append using "aggregated_bf_itn.dta"
*Remove variables we don't need anymore
drop sum1 sum2 prop_fed or1 or2

/* B. At the moment, 'itn' is a string variable. We need a 
factor variable for the regression model */
encode(itn), generate(itn2)
replace itn2=. if day==.
levelsof(itn2)
*See how the itn2 levels correspond to itn*
label list itn2

/* C. For the regression model we use the function 'blogit', 
which fits a logistic regression model to aggregated data */
blogit tot_bf total i.itn2 i.hut i.sleeper i.day i.wash

/* Calculate the odds ratio (OR) and 95% CI for the candidate net,
 which is level 2 of itn2 */
gen or_model = exp(_b[_outcome:2.itn2])
gen or_model_lower = exp(_b[_outcome:2.itn2] - 1.96* _se[_outcome:2.itn2])
gen or_model_upper = exp(_b[_outcome:2.itn2] + 1.96* _se[_outcome:2.itn2])

*Alternatively, we could have asked Stata to calculate the ORs for us, like this:*
blogit tot_bf total i.itn2 i.hut i.sleeper i.day i.wash, or

*Recall the non-inferiority margin. Use this to compare with the OR calculated above
list if itn=="Active_comparator" & missing(day)

/* D. Now check whether the candidate net is superior to the
standard comparator (level 4 of 'itn2').  */
blogit tot_bf total ib4.itn2 i.hut i.sleeper i.day i.wash

/* ############################################################## 
							IRS Example     
############################################################## */


/* As the procedure for IRS is extremely similar to that carried out for ITNs,
We won't repeat all the details here. But we will show a quick example. 
Using the dataset loaded below, we'll test whether a candidate IRS product is 
non-inferior to an active comparator IRS product, in terms of mosquito mortality.
Here, we'll compare how the products performed when sprayed on a mud wall. 
In the dataset, products were assessed on both mud and concrete. An untreated
was also included */

clear
import delimited "example_dataset_IRS.csv"

tab tot_dead
summarize tot_dead
levelsof(hut)
levelsof(sleeper)
levelsof(treatment)

/* Determine the non-inferiority margin (NIM)
NIM: mosquito mortality induced by the candidate product should be no more than 7% less
 than that induced by the active comparator product.
First, calculate the (unadjusted) mosquito mortality in each trial arm */
collapse (sum) sum1=tot_dead sum2=total, by(treatment)
gen prop_dead = sum1/sum2
list
gen or1 = (prop_dead - 0.07)/(1-prop_dead + 0.07)
gen or2 = (prop_dead)/(1-prop_dead)
*Calculate the odds-ratio (OR) for the NIM
gen nim = or1/or2
list
save "aggregated_mortality_IRS.dta"

/* ############################################################## 
		Mosquito Mortality (products applied on mud walls)   
 ############################################################## */

/* A. Reload dataset to carry out the regression*/
clear
import delimited "example_dataset_IRS.csv"
append using "aggregated_mortality_IRS.dta"
*Remove variables we don't need anymore
drop sum1 sum2 prop_dead or1 or2

/* B. At the moment, 'treatment' is a string variable. We need a 
factor variable for the regression model */
encode(treatment), generate(treatment2)
replace treatment2=. if day==.
levelsof(treatment2)
*See how the treatment2 levels correspond to treatment*
label list treatment2

/* C. For the regression model we use the function 'blogit', 
which fits a logistic regression model to aggregated data.
As treatments are not rotated around huts, 
the variable 'hut' is less informative here.
We want the 2nd level of treatment2 as baseline (active comparator on mud) */
blogit tot_dead total ib2.treatment2 i.sleeper i.day

*Here is how the model is stored in Stata's memory*
ereturn list

*Calculate the odds ratio (OR) and 95% CI for the unwashed candidate product
gen or_model = exp(_b[_outcome:4.treatment2])
gen or_model_lower = exp(_b[_outcome:4.treatment2] - 1.96* _se[_outcome:4.treatment2])
gen or_model_upper = exp(_b[_outcome:4.treatment2] + 1.96* _se[_outcome:4.treatment2])

*Alternatively, we could have asked Stata to calculate the ORs for us, like this:*
blogit tot_dead total ib2.treatment2 i.sleeper i.day, or

*Recall the non-inferiority margin. Use this to compare with the OR calculated above
list if treatment=="Active_comp_mud" & missing(day)


/* End of IRS example*/


/* The code in this section shows how to visualise the non-inferiority
assessment in STATA, in a similar way to that shown in the R tutorial. 
Thanks to John Bradley for his assistance with this */

clear
set scheme s1color

*plot showing NI margins
set obs 800
gen mm = (_n/10) + 19.9 
gen NI = ((mm-7)/(107-mm)) / (mm/(100-mm))  

twoway (area mm NI, color(red*0.4)), ///
ytitle("Mortality for active comparator (%)") xtitle("Odds ratio") /// 
ylabel(20(10)100, angle(horizontal)) xlabel(0.1(.1)2.0) ///
plotregion(m(b=0 l=0 r=0 t=0) )  plotregion(color(green*0.4)) ///
xline(1, lpatter(dash) lcolor(red*0.6))

*with exmaple OR and CI plotted
gen OR = 0.95 
gen OR_lb = 0.8
gen OR_ub = 1.13
gen cont_mort = 60 

twoway (area mm NI, color(red*0.4)) ///
(scatter cont_mort OR, color(black)) ///
(rcap OR_ub OR_lb cont_mort , horizontal color(black)), ///
ytitle("Mortality for active comparator (%)") xtitle("Odds ratio") /// 
ylabel(20(10)100, angle(horizontal)) xlabel(0.1(.1)2.0) ///
plotregion(m(b=0 l=0 r=0 t=0) )  plotregion(color(green*0.4)) ///
legend(off) xline(1, lpatter(dash) lcolor(red*0.6))


*plotting on a reasonable scale for this OR - this takes a little playing around with to get right...
replace NI = 0.5 if mm > 79.99 & mm < 80.01
set obs 801
replace mm = 79.999 in 801
replace NI = ((79.999-7)/(107-79.999)) / (79.999/(100-79.999)) in 801
sort mm
twoway (area mm NI, color(red*0.4)) ///
(scatter cont_mort OR, color(black)) ///
(rcap OR_ub OR_lb cont_mort , horizontal color(black)) ///
if mm <= 80 & mm >= 30, ///
ytitle("Mortality for active comparator (%)") xtitle("Odds ratio") /// 
ylabel(30(10)80, angle(horizontal)) xlabel(0.5(.1)1.3) ///
plotregion(m(b=0 l=0 r=0 t=0) )  plotregion(color(green*0.4)) ///
legend(off) xline(1, lpatter(dash) lcolor(red*0.6))





