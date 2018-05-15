clear all

//Set various loop parameters etc:

set seed 123
*Number of simulations
local S = 1000

//Number of obs per simulation
local n_values 100 , 500, 1000 

//Note relative power is defined below because it would be a PITA to define it up here like I should


//Define a program that runs a single iteration of the simulation
cap program drop simulate
program  define simulate , rclass
	syntax, n(real)  beta_0(real)
		clear

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
		lincom /beta
		matrix b = b, `r(se)' 
		
		test /beta = 0
		matrix b = b, `r(p)'

		return matrix b b
end





//Loop over simulation specifications:
//set outcomes table matrix
matrix outcomes_tab = J(1, 7,.)
foreach n of numlist `n_values' { 
	local relative_power = .25/`n'
	foreach beta_0 of numlist `relative_power', .25, .10 {
		//Loop over the program to simulate the data as requested
		matrix  outcomes = J(1,4,.)
		forvalues s = 1/`S' {
			simulate, n(`n') beta_0(`beta_0')
			//matrix b = r(b)
			matrix outcomes = outcomes \ r(b)
		}

		//turn simulation into data to examine
		clear
		svmat outcomes
		//drop the first row that initialized the data
		keep if !missing(outcomes1)

		//Bias_Bar
		rename (outcomes1 outcomes2 outcomes3 outcomes4) (beta_hat sigma_Hat beta_SE p)
		sum beta_hat
		local beta_hat_bar = `r(mean)'
		local bias_bar = `beta_hat_bar' - `beta_0'

		//SSD
		generate u_hat_square = (beta_hat - `beta_hat_bar')^2
		sum u_hat_square 
		local SSD = `r(mean)'

		//SE/SSD
		sum beta_SE
		local SE_bar = `r(mean)'
		local SE_SSD = `SE_bar'/`SSD'

		//Test Size
		gen reject_05 = p < .05
		gen reject_10 = p < .10

		sum reject_05
		local p_bar_05 = `r(mean)'

		sum reject_10
		local p_bar_10 = `r(mean)'

		//Combine into a single row, and append to outcomes_tab
		matrix line = (`n', `beta_0', `bias_bar', `SSD', `SE_SSD', `p_bar_05', `p_bar_10')	
		matrix outcomes_tab = outcomes_tab \ line
	}
}
matrix colnames outcomes_tab = N  beta_0 bias_bar SSD SE_SSD p_bar_05 p_bar_10

 



