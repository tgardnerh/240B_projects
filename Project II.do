log using bcuse jtrain, clear
regress lscrap d88 d89 grant grant_1
xtset fcode
xtreg lscrap d88 d89 grant grant_1, re
estimates store regls

xtreg lscrap d88 d89 grant grant_1, fe
estimates store fe

capture xtreg lscrap d88 d89 grant grant_1 union, fe
display "`_rc'"

#delimit ;
display 
"Calling the estimates from the FE model $\hat{\beta_0}}$ 
and the results from the RE model $\hat{\beta_1}}$,
 the Hausman statisic to assess whether RE is a reasonable assumption is 
 
$$(\hat{\beta_1}}-\hat{\beta_0}}'(\hat{\text{Var}(\hat{\beta_1}})}+\hat{\text{Var}(\hat{\beta_0}})})(\hat{\beta_1}}-\hat{\beta_0}})$$

";

#delimit cr

hausman fe regls


xtreg lscrap d88 d89 grant grant_1, re vce(robust)
estimates store regls_robust

xtreg lscrap d88 d89 grant grant_1, fe vce(robust)
estimates store fe_robust

#delimit ;
display 
"Standard Errors of the robust estimates are significantly larger, suggesting that
the GLS model is efficient, if it is consistent."
;
#delimit cr



