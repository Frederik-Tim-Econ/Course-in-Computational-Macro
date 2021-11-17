# import relevant packages
import numpy as np
import Tools as tools

def solve(par): # Solves the model using backwards induction and EGM
    # Preallocating
    Cstar = np.zeros([par['T'],par['gridsize_a']+1])
    𝒢_w   = np.zeros([par['T'],par['gridsize_a']+1])    
    
    # 1) Loop over time
    for t in range(par['T']-1,-1,-1): 
        print(t)       
        if t == par['T'] - 1:
            # Use an exogenous grid in the last period
            𝒢_w[t,:]   = par['𝒢_w_det']
            Cstar[t,:] = np.array(𝒢_w[t,:]).T              
        else:            
            # 2) Loop over savings in exogenous grid
            for ia,a in enumerate(par['𝒢_a']):
                # Cash-on-hand tomorrow
                w_plus = par['l'][t+1]*par['Y'] + par['R']*a
                
                # Find consumption tomorrow
                C_plus = tools.interp_linear_1d(𝒢_w[t+1,:],Cstar[t+1,:],w_plus)
                
                # Stochastic Euler
                Cstar[t,ia+1] = (par['R']*par['β']*(par['ω'] @ C_plus**(-par['ρ'])))**(-1/par['ρ'])
                
                # Updating the Endogenous Grid
                𝒢_w[t,ia+1] = a + Cstar[t,ia+1]
    return 𝒢_w,Cstar