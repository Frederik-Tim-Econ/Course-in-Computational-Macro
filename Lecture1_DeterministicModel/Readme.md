# Lecture 1: A Deterministic General-Equilibrium Life-cycle Model

We consider a deterministic finite-horizon model where individuals live for 𝑇 periods and work for 𝑟 periods.
We first solve for optimal consumption and savings closed-form for given factor prices.
We then use numerical solution methods to determine a steady-state where factor prices are consistent with the level of savings and, by extension, capital.

## OVERVIEW OF FILES
 * Lecutre_1.ipynb:
	* A notebook containing the description of the model and its implementation
	* Defines the exogenous parameter values in dictionary "par".
	* Defines and calls the functions: 
		* solve(): 	Solves for individual decisions for given factor prices
		* objective():	Defines a loss-function to be minimized over a criterion for the interest rate
 * Exercise.ipynb:
	* A notebook stating the exercise
	* You are asked to introduce a pension scheme to the model
 * Exercise_Solution.ipynb:
	* A notebook with the solution code to the exercise
	* Defines and calls the functions: 
		* solve_paygo(): 	Solves for individual decisions for given factor prices in a model with a pension scheme
		* objective_paygo():	Defines a loss-function to be minimized over a criterion for the interest rate
 