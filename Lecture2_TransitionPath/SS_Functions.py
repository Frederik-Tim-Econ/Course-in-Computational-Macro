import numpy as np
import scipy.optimize as optimize


# Solves for steady state capital given pension contribution rate 
def Kss(τ,par):
    
    # Unpacking Parameters
    α = par['α']
    δ = par['δ']
    L = par['L']
    
    r_lb = 1e-4 - δ # lower bound r
    K_ub = L*((r_lb+δ)/α)**(1/(α-1)) # upper bound K(r_lb)
    sol = optimize.minimize_scalar(objective,bounds=[0,K_ub],args=(τ,par),method='bounded',options={'xatol': 1e-4, 'maxiter': 10000})
    K_ss = sol.x
    
    return K_ss

# new objective function calling solve_paygo
def objective(K_guess,τ,par): # Solves the SS
    
    C,A,K_implied = solve(K_guess,τ,par) 
    
    # STEP 8: Check distance between K_guess - K_implied (Note we define the loss as the squared difference)
    loss = (K_guess - K_implied)**2  
    return loss


def solve(K_guess,τ,par): # Solve for consumption, savings and aggregate capital given a guess for r
    
    
    # Unpacking Parameters
    T = par['T']
    R = par['R']
    α = par['α']
    ρ = par['ρ']
    δ = par['δ']
    β = par['β']
    l = par['l']
    L = par['L']
    
    
   
    # STEP 2: Solve for wage given guessed interest rate
    r = α*(K_guess/L)**(α-1) - δ
    w = (1 - α)*(K_guess/L)**α
    
    

    # STEP 3: Solve for first periond consumption (we have a closed form in (1))
    b = τ*w*sum(l) / (T-R) # solve for the steady state pension benefit
    C1 = sum(( (1-τ)*l*w +  (1-l)*b ) / (1+r)**np.linspace(1,T,T)) /sum((β*(1+r))**((np.linspace(1,T,T)-1)/ρ)/(1+r)**np.linspace(1,T,T))
    
    # STEP 4: Solve for the whole consumption path using the Euler equation
    C = C1*(β*(1+r))**((np.arange(T))/ρ)

    # STEP 5: Solve for the whole saviongs path using the budget constraint
    A = np.zeros(T) # preallocate storage
    A[0] = (1-τ)*w*l[0] - C[0] # solve for first period savings given no initial wealth
    for t in range(T): # solve the whole savings path
        if t>0:
            A[t] = (1-τ)*w*l[t] + (1-l[t])*b + (1 + r)*A[t-1] - C[t]

    # STEP 5: Compute implied aggregate capital by summing over savings path
    K = sum(A)
     
    return C,A,K

