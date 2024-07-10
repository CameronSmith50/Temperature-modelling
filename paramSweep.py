# Functions which will enable the generation and plotting of datasets with
# random parameter sets to determine patterns in when resistance is selected
# for
#
# Author: Cameron Smith

#%% Packages and files

from ecoSimND import nondimensionalSim
from Initialisation import initNonDim
import Lookup as LU
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import sys
import os

#%% Function to generate data

def generateData(pars, n=100, save=False):
    '''
    Function which will generate a .csv file which will contain the random
    parameter set, together with the required outputs such as the state
    variables at steady state.
    The main input is pars, which should come from the function initNonDim. Note
    that several of these parameters will be overwritten in the generation of
    the dataset.
    n is the number of simulations to run and save determines if the dataset
    should be saved (True) or output only (False)
    '''

    # Decide on ranges for the various parameters
    rRange = (1.0, 5.0) # A transformed variable for r. This will be the 
                        # change, with an extra parameter deciding if the 
                        # change is smaller or larger
    phiWRange = (0.0, 1.0) # Wildtype fecundity loss
    phiMRange = (0.0, 1.0) # Mutant fecundity loss 
    alphaWRange = (0.0, 0.5) # Wildtype virulence
    gammaRange = (0.0, 20.0) # Dilution effects

    # Create an empty dataframe
    cols = ['r', 'phiW', 'phiM', 'alphaW', 'gamma', 'uW', 'vW', 
                'uM', 'vM', 'n', 'P']
    df = pd.DataFrame(columns = cols)

    # Draw the random variables and transform if required
    rVec = [n if np.random.rand() < 0.5 else 1/n for n in \
            ((rRange[1]-rRange[0])*np.random.rand(n) + rRange[0])]
    phiWVec = (phiWRange[1]-phiWRange[0])*np.random.rand(n) + phiWRange[0]
    phiMVec = (phiMRange[1]-phiMRange[0])*np.random.rand(n) + phiMRange[0]
    alphaWVec = (alphaWRange[1]-alphaWRange[0])*np.random.rand(n) + \
            alphaWRange[0]
    gammaVec = (gammaRange[1]-gammaRange[0])*np.random.rand(n) + \
            gammaRange[0]

    # Loop through the n variable
    for kk in tqdm(range(n), leave = False, ncols = 50):
        
        # Create an empty dictionary for this repeat and populate with the
        # random parameters
        dfDict = {}
        dfDict['r'] = [rVec[kk]]
        dfDict['phiW'] = [phiWVec[kk]]
        dfDict['phiM'] = [phiMVec[kk]]
        dfDict['alphaW'] = [alphaWVec[kk]]
        dfDict['gamma'] = [gammaVec[kk]]

        # Update pars
        pars['growthRatio'] = rVec[kk]
        pars['fecunW'] = phiWVec[kk]
        pars['fecunM'] = phiMVec[kk]
        pars['virNDW'] = alphaWVec[kk]
        pars['pathDilND'] = gammaVec[kk]

        # Run a simulation
        sim = nondimensionalSim(pars)
        sim.solveSystem()

        # Add to the dictionary with the outputs
        dfDict['uW'] = [sim.SS[0] if np.abs(sim.SS[0])>1e-3 else 0.0] 
        dfDict['vW'] = [sim.SS[1] if np.abs(sim.SS[1])>1e-3 else 0.0] 
        dfDict['uM'] = [sim.SS[2] if np.abs(sim.SS[2])>1e-3 else 0.0] 
        dfDict['vM'] = [sim.SS[3] if np.abs(sim.SS[3])>1e-3 else 0.0] 
        dfDict['n'] = [np.sum(sim.SS[:4]) if np.abs(np.sum(sim.SS[:4]))>1e-3 \
                    else 0.0]
        dfDict['P'] = [sim.SS[4] if np.abs(sim.SS[4])>1e-3 else 0.0] 

        # Append this repeat to the wider dataFrame
        df1 = pd.DataFrame().from_dict(dfDict)
        df = pd.concat([df, df1], ignore_index=True)

    # Save if necessary
    if save:

        # Generate a save name
        saveName = LU.createSaveName(f'SCAT')
        
        # Create the two directories in ./Data and ./Figures
        os.mkdir(f'./Data/{saveName}/')
        
        # Create a parameters file and store it in the Data directory
        del pars['growthRatio']
        del pars['fecunW']
        del pars['fecunM']
        del pars['virNDW']
        del pars['pathDilND']

        with open(f'./Data/{saveName}/' + 'parameters.txt', 'w') as file:
            file.write('Parameters:\n')
            file.write('----------\n')
            file.write('\n')
            for key, val in pars.items():
                file.write(f'{key: <13}: {val:.3f}\n')
        
        # Save the dataframe
        df.to_csv(f'./Data/{saveName}/dataset.csv', index=False)

        # Add an entry to the lookup table
        LU.addEntrytoLookup(saveName, pars, 'SCAT')

    else:
    
        return(df)
    
#%% Plotting

def plotDataTest(folderName, save = False):
    '''
    Code to generate a figure of the scatter plots. folderName should be the
    folder name containing data. If save is False, we will just plot the figure.
    If True, the option to overwrite the default savename or create a new figure
    will be given.
    '''

    # Firstly, check if the folderName exists
    if not os.path.isdir(f'./Data/{folderName}'):
        print('Folder does not exist in the Data directory.')
        return
    
    # Load the dataset
    df = pd.read_csv(f'./Data/{folderName}/dataset.csv')

    # Generate some new columns in the dataframe from old columns
    df['xW'] = (df['uW'] + df['vW']) # Density of wildtype
    df['xM'] = (df['uM'] + df['vM']) # Density of mutant
    df['chi'] = (df['xW'] - df['xM'])/df['n'] # Ratio of wildtype to mutant

    # Generate a figure
    figWidth = 20
    figHeight = 20
    fsize = 20
    fig = plt.figure(figsize=(figWidth, figHeight))
    axr = fig.add_subplot(321)
    axFecun = fig.add_subplot(322)
    axDil = fig.add_subplot(323)
    axChi = fig.add_subplot(324)
    axComp = fig.add_subplot(325)
    axTemp = fig.add_subplot(326)

    # Plot 1: Growth ratio vs chi
    axr.scatter(np.log10(df['r']), df['chi'], c=df['chi'], 
                    s=5, cmap='coolwarm_r', vmin=-1.0, vmax=1.0)
    axr.set_xlabel(r'$log_{10}(r)$', fontsize=fsize)
    axr.set_xticks([np.log10(min(df['r'])), 0.0, np.log10(max(df['r']))])
    axr.set_xticklabels([f"{np.log10(min(df['r'])):.2f}", '0', 
                         f"{np.log10(max(df['r'])):.2f}"], fontsize=fsize)
    axr.set_ylabel(r'$\chi$', fontsize=fsize)
    axr.set_yticks([-1, 0, 1])
    axr.set_yticklabels(['-1', '0', '1'], fontsize=fsize)

    # Plot 2: Fecundity losses vs chi
    axFecun.scatter(df['phiW'], df['phiM'], c=df['chi'], s=5, cmap='coolwarm_r',
                    vmin=-1.0, vmax=1.0)
    axFecun.set_xlabel(r'Wildtype fecundity, $\varphi_W$', fontsize=fsize)
    axFecun.set_xticks([0, 0.5, 1])
    axFecun.set_xticklabels(['0', '0.5', '1'], fontsize=fsize)
    axFecun.set_ylabel(r'Mutant fecundity, $\varphi_M$', fontsize=fsize)
    axFecun.set_yticks([0, 0.5, 1])
    axFecun.set_yticklabels(['0', '0.5', '1'], fontsize=fsize)

    # Plot 3: Fecundity ratio and growth ratio vs chi
    axDil.scatter(df['gamma'], df['chi'], c=df['chi'], s=5, cmap='coolwarm_r',
                    vmin=-1.0, vmax=1.0)
    axDil.set_xlabel(r'Dilution effects, $\gamma$', fontsize=fsize)
    axDil.set_xticks([min(df['gamma']), 
                      0.5*(min(df['gamma']) + max(df['gamma'])), 
                      max(df['gamma'])])
    axDil.set_xticklabels([f"{min(df['gamma']):.1f}", 
                      f"{0.5*(min(df['gamma']) + max(df['gamma'])):.1f}", 
                      f"{max(df['gamma']):.1f}"], fontsize=fsize)
    axDil.set_ylabel(r'$\chi$', fontsize=fsize)
    axDil.set_yticks([-1, 0, 1])
    axDil.set_yticklabels(['-1', '0', '1'], fontsize=fsize)
    
    # Plot 4: Histogram for the chi parameter
    my_cmap = plt.get_cmap("coolwarm_r")
    rescale = lambda y: (y - np.min(y)) / (np.max(y) - np.min(y))
    bins = np.linspace(-1, 1, 21)
    histCounts = np.histogram(df['chi'], bins)[0]/len(df)
    binmps = 0.5*(bins[:-1]+bins[1:])
    axChi.bar(binmps, histCounts, width=bins[1]-bins[0],
                    color=my_cmap(rescale(binmps)))
    
    axChi.set_xlabel(r'$\chi$', fontsize=fsize)
    axChi.set_xticks([-1, 0, 1])
    axChi.set_xticklabels(['-1', '0', '1'], fontsize=fsize)

    # Plot 5: Different bars depending on the dilution effects
    maxGamma = max(df['gamma'])
    dfLowGamma = df[df['gamma'] < maxGamma*0.1]
    dfHighGamma = df[df['gamma'] > maxGamma*0.9]
    bins = np.linspace(-1, 1, 11)
    lowHist = np.histogram(dfLowGamma['chi'], bins)[0]/len(dfLowGamma + \
                        dfHighGamma)
    highHist = np.histogram(dfHighGamma['chi'], bins)[0]/len(dfLowGamma + \
                        dfHighGamma)
    binmps = (bins[:-1] + bins[1:])*0.5
    db = bins[1]-bins[0]
    axComp.bar(binmps, lowHist, width=db/2, color='k', align='edge',
               label = 'low dilution')
    axComp.bar(binmps, highHist, width=-db/2, color='g', align='edge',
               label = 'High dilution')
    axComp.legend()



    axTemp.scatter(df['n'], df['chi'], c=df['chi'], s=5, cmap='coolwarm_r',
                   vmin=-1, vmax=1)

    # Save if necessary
    if save:
        if not os.path.isdir(f'./Figures/{folderName}'):
            os.mkdir(f'./Figures/{folderName}')
            fig.savefig(f'./Figures/{folderName}/results.png', 
                    bbox_inches='tight')
        else:
            newFig = input('Would you like a new plot. If no, this will '
                    'overwrite the "scatter*" plots already in the folder.'
                     ' [y/n]: ')
            if newFig == 'y':
                newName = input('What name would you like (no extension)? ')
                fig.savefig(f'./Figures/{folderName}/{newName}.png', 
                    bbox_inches='tight')
            elif newFig == 'n':
                fig.savefig(f'./Figures/{folderName}/results.png', 
                    bbox_inches='tight')
            else:
                print('Invalid response, exiting code.')
                return
    else:
        plt.show()

def plotData(folderName, save = False):
    '''
    Code to generate a figure of the scatter plots. folderName should be the
    folder name containing data. If save is False, we will just plot the figure.
    If True, the option to overwrite the default savename or create a new figure
    will be given.
    '''

    # Firstly, check if the folderName exists
    if not os.path.isdir(f'./Data/{folderName}'):
        print('Folder does not exist in the Data directory.')
        return
    
    # Load the dataset
    df = pd.read_csv(f'./Data/{folderName}/dataset.csv')

    # Generate some new columns in the dataframe from old columns
    df['xW'] = (df['uW'] + df['vW']) # Density of wildtype
    df['xM'] = (df['uM'] + df['vM']) # Density of mutant
    df['chi'] = (df['xW'] - df['xM'])/df['n'] # Ratio of wildtype to mutant

    # Generate a figure
    figWidth = 20
    figHeight = 20
    fsize = 20
    fig = plt.figure(figsize=(figWidth, figHeight))
    axr = fig.add_subplot(221)
    axDil = fig.add_subplot(222)
    axChi = fig.add_subplot(223)
    axComp = fig.add_subplot(224)

    # Plot 1: Growth ratio vs chi
    axr.scatter(np.log10(df['r']), df['chi'], c=df['chi'], 
                    s=20, cmap='coolwarm_r', vmin=-1.0, vmax=1.0)
    axr.set_xlim([np.log10(min(df['r'])), np.log10(max(df['r']))])
    axr.set_xlabel(r'$log_{10}(\lambda)$', fontsize=fsize)
    axr.set_xticks([np.log10(min(df['r'])), 0.0, np.log10(max(df['r']))])
    axr.set_xticklabels([f"{np.log10(min(df['r'])):.2f}", '0', 
                         f"{np.log10(max(df['r'])):.2f}"], fontsize=fsize)
    axr.set_ylim([-1, 1])
    axr.set_ylabel(r'Susceptible genotype ratio, $\chi$', fontsize=fsize)
    axr.set_yticks([-1, 0, 1])
    axr.set_yticklabels(['-1', '0', '1'], fontsize=fsize)
    axr.text(-0.7, 1.05, 'A.', fontsize=fsize)

    # Plot 2: Fecundity ratio and growth ratio vs chi
    axDil.scatter(df['gamma'], df['chi'], c=df['chi'], s=20, cmap='coolwarm_r',
                    vmin=-1.0, vmax=1.0)
    
    maxGamma = max(df['gamma'])
    axDil.set_xlim([0, maxGamma])
    axDil.annotate('', xy=(0.95, -0.025), xycoords='axes fraction', 
                   xytext=(0.05, -0.025), 
                   arrowprops=dict(color='k'))
    axDil.text(10.0, -1.125, "Increasing dilution effects", fontsize=fsize,
               ha='center', va='center')
    axDil.set_xticks([])
    axDil.set_xticklabels([])
    axDil.set_ylim([-1, 1])
    axDil.set_ylabel(r'Susceptible genotype ratio, $\chi$', fontsize=fsize)
    axDil.set_yticks([-1, 0, 1])
    axDil.set_yticklabels(['-1', '0', '1'], fontsize=fsize)
    axDil.text(0.0, 1.05, 'B.', fontsize=fsize)
    
    # Plot 3: Histogram for the chi parameter
    my_cmap = plt.get_cmap("coolwarm_r")
    rescale = lambda y: (y - np.min(y)) / (np.max(y) - np.min(y))
    bins = np.linspace(-1, 1, 11)
    histCounts = np.histogram(df['chi'], bins)[0]/len(df)
    binmps = 0.5*(bins[:-1]+bins[1:])
    axChi.barh(binmps, histCounts, bins[1]-bins[0], edgecolor='k',
                    color=my_cmap(rescale(binmps)))
    maxx = np.ceil(np.max(histCounts)*20)/20
    
    axChi.set_xticks([0.0, maxx/2, maxx])
    axChi.set_xticklabels(['0', f'{maxx/2:.2f}', f'{maxx:.2f}'], fontsize=fsize)
    axChi.set_ylim([-1, 1])
    axChi.set_ylabel(r'Susceptible genotype ratio, $\chi$', fontsize=fsize)
    axChi.set_yticks([-1, 0, 1])
    axChi.set_yticklabels(['-1', '0', '1'], fontsize=fsize)
    axChi.text(0.0, 1.05, 'C.', fontsize=fsize)

    # Plot 4: Different bars depending on the dilution effects
    dfLowGamma = df[df['gamma'] < maxGamma*0.1]
    dfHighGamma = df[df['gamma'] > maxGamma*0.9]
    lowHist = np.histogram(dfLowGamma['chi'], bins)[0]/len(df)
    highHist = np.histogram(dfHighGamma['chi'], bins)[0]/len(df)
    binmps = (bins[:-1] + bins[1:])*0.5
    db = bins[1]-bins[0]
    axComp.barh(binmps, lowHist, db/2, color=my_cmap(rescale(binmps)), 
                edgecolor='k', align='edge', hatch='x')
    axComp.barh(binmps, highHist, -db/2, color=my_cmap(rescale(binmps)), 
                edgecolor='k', align='edge')
    axComp.barh(-5, 1e-3, 1e-3, color='w', edgecolor='k', hatch='x', 
                label='Low dilution')
    axComp.barh(-5, 1e-3, 1e-3, color='w', edgecolor='k',
                label='High dilution')
    maxx = np.ceil(np.max([np.max(highHist), np.max(lowHist)])*100)/100
    
    axComp.set_xticks([0.0, maxx/2, maxx])
    axComp.set_xticklabels(['0', f'{maxx/2:.2f}', f'{maxx:.2f}'],
                fontsize=fsize)
    axComp.set_ylim([-1, 1])
    axComp.set_ylabel(r'Susceptible genotype ratio, $\chi$', fontsize=fsize)
    axComp.set_yticks([-1, 0, 1])
    axComp.set_yticklabels(['-1', '0', '1'], fontsize=fsize)
    axComp.text(0.0, 1.05, 'D.', fontsize=fsize)
    axComp.legend(fontsize=fsize)

    # Save if necessary
    if save:
        if not os.path.isdir(f'./Figures/{folderName}'):
            os.mkdir(f'./Figures/{folderName}')
            fig.savefig(f'./Figures/{folderName}/results.png', 
                    bbox_inches='tight')
        else:
            newFig = input('Would you like a new plot. If no, this will '
                    'overwrite the "scatter*" plots already in the folder.'
                     ' [y/n]: ')
            if newFig == 'y':
                newName = input('What name would you like (no extension)? ')
                fig.savefig(f'./Figures/{folderName}/{newName}.png', 
                    bbox_inches='tight')
            elif newFig == 'n':
                fig.savefig(f'./Figures/{folderName}/results.png', 
                    bbox_inches='tight')
            else:
                print('Invalid response, exiting code.')
                return
    else:
        plt.show()

if __name__ == '__main__':

    # Testing
    if len(sys.argv) < 2:
        pars = initNonDim()
        print(generateData(pars, n=100))
    
    # Saving data
    elif sys.argv[1] == 'save':
        pars = initNonDim(K=1,
                          tEco=2000)
        generateData(pars, n = 10000, save = True)

    # Plotting in window
    elif sys.argv[1] == 'plot':
        plotData(sys.argv[2])

    # Save plotting
    elif sys.argv[1] == 'saveplot':
        plotData(sys.argv[2], save=True)
