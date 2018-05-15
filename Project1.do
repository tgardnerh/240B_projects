clear all

program define censored_model
	args lnf beta sigma
	#delimit ;
	quietly replace `lnf'=censored_d*(ln(normalden(($ML y-`beta'*x)/`sigma'))- ln(`sigma')) 
		+ (1-censored_d)*(1-ln(normal(($ML y-`beta'*x)/`sigma'))- ln(`sigma')) 
		;
	#delimit cr
end

program define normal
args lnf mu sigma
quietly replace `lnf' = ln(normal(($ML_y1-`mu')/`sigma') )-ln(`sigma')
end

program define bern
args lnf mu 
quietly replace `lnf' = $ML_y1*ln(`mu') + (1-$ML_y1)*ln(`mu')
end


cap program drop simulate

program  define simulate , rclass

	syntax, n(real)  beta_0(real)
		

		set obs `n'
		generate x = runiform(-5,5)
		generate u = rnormal(0,1)
		generate y_star = `beta_0'*x + u

		generate y = min(y_star, 1)
		generate censored_d = y == 1
	
	
		#delimit ;
		mlexp (
			censored_d*( lnnormalden((y-{beta}*x)/{sigma})- ln({sigma}) )
			+ (1-censored_d)*ln(1-normal((y-{beta}*x)/{sigma}))- ln({sigma})
			) 
			;
		#delimit cr
		matrix define b = e(b)
		return matrix b b

end

simulate, n(1000) beta_0(0)






