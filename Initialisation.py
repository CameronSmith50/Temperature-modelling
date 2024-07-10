# Code which defines the parameters for either a dimensional or non-dimensional
# system.
# 
# Author: Cameron Smith
# Date created: 09/04/2024

#%% Dimensional system

def initDim(
        rW = 2.5, # Dimension: per unit time
        rM = 3.7, # Dimension: per unit time
        phiW = 0.73, # Dimension: dimensionless
        phiM = 0.19, # Dimension: dimensionless
        beta = 0.5, # Dimension: per pathogen per unit time
        theta = 0.1, # Dimension: per unit time
        sigmaW = 0.13, # Dimension: per unit time
        sigmaM = 0.0, # Dimension: per unit time
        q = 1000.0, # Dimension: number of hosts
        a = 1e-3, # Dimension: pathgogens per infected individual per unit time
        rho = 1.0, # Dimension: per unit time
        K = 1e3, # Dimension: number of pathogens
        delta = 0.01, # Dimension: per unit time
        eta = 1e-4, # Dimension: pathogens per susceptible individual
        gamma = 0.0, # Dimension: per resistant individual per time
        tEco = 100.0 # Dimension: unit time
    ):
    '''
    Initial parameters for the fully dimensional system.
    '''

    # Initialise the dictionary
    pars = {}

    # Growth rates
    pars['growthW'] = rW
    pars['growthM'] = rM

    # Loss in fecundity while infected
    pars['fecunW'] = phiW
    pars['fecunM'] = phiM

    # Transmission
    pars['transmission'] = beta

    # Natural mortality
    pars['natMort'] = theta

    # Virulences
    pars['virW'] = sigmaW
    pars['virM'] = sigmaM

    # Carrying capacity
    pars['comp'] = q

    # Pathogen parameters
    pars['pathShed'] = a
    pars['pathDecay'] = delta
    pars['pathLoss'] = eta
    pars['pathDil'] = gamma
    pars['pathGrow'] = rho
    pars['pathCap'] = K

    # Final time of ecological simulation
    pars['tEco'] = tEco

    return(pars)

#%% Non-dimensional system

def initNonDim(
        r = 0.1,
        phiW = 0.7,
        phiM = 0.2,
        d = 0.1,
        alphaW = 0.5,
        alphaM = 0.0,
        a = 1e-3,
        rho = 1.0,
        K = 1,
        delta = 0.01,
        eta = 1e-4,
        gamma = 0.0,
        tEco = 200.0
    ):
    '''
    Initial parameters for the non-dimensional system
    '''

    # Initialise the dictionary
    pars = {}

    # Ratio of growth rates
    pars['growthRatio'] = r

    # Loss of fecundity on infection
    pars['fecunW'] = phiW
    pars['fecunM'] = phiM

    # ND natural mortality
    pars['natMortND'] = d

    # ND virulences
    pars['virNDW'] = alphaW
    pars['virNDM'] = alphaM

    # ND pathogen parameters
    pars['pathShedND'] = a
    pars['pathDecayND'] = delta
    pars['pathLossND'] = eta
    pars['pathDilND'] = gamma
    pars['pathCapND'] = K
    pars['pathGrowND'] = rho

    # Final time of ecological simulation
    pars['tEco'] = tEco

    return(pars)

# Output the initial parameters if the file is run in the terminal
if __name__ == '__main__':

    # Dimensional
    print('')
    print('Dimensional system')
    print('------------------')
    parsDim = initDim()
    for key, val in parsDim.items():
        print(f'{key}: {val:.3f}')

    # Non-dimensional
    print('')
    print('Non-dimensional system')
    print('----------------------')
    parsND = initNonDim()
    for key, val in parsND.items():
        print(f'{key}: {val:.3f}')

    
