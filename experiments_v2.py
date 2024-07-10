# Functions which will enable the generation and plotting of datasets with
# random parameter sets which more closely align to the experiments
#
# Author: Cameron Smith

#%% Packages and files

from ecoSimDim import dimensionalSim
from Initialisation import initDim
import Lookup as LU
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import sys
import os
import scipy.stats as sp

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

    # Create an empty dataframe
    cols = ['Experiment', 'rW', 'rM', 'phiW', 'phiM', 'sigmaW', 'gamma', 'SW',
            'IW', 'SM', 'IM', 'N', 'P']
    df = pd.DataFrame(columns = cols)

    # Draw the random variables and transform if required
    # Create several empty np arrays
    rWVec = np.zeros(4*(n//4))
    rMVec = np.zeros(4*(n//4))
    phiWVec = np.zeros(4*(n//4))
    phiMVec = np.zeros(4*(n//4))
    dRatio = np.zeros(4*(n//4))
    sigmaWVec = np.zeros(4*(n//4))

    # Experiment 20-20
    for kk in range(n//4):
        rWVec[kk] = np.random.normal(52.3/24, 23.3/576)
        # rWVec[kk] = np.random.rand()*4 + 0.5
        rMVec[kk] = np.random.normal(77/24, 18.9/576)
        while phiWVec[kk] <= 0 or phiWVec[kk] >= 1:
            phiWVec[kk] = np.random.normal(0.728, 0.233)
            # phiWVec[kk] = np.random.rand()
        while phiMVec[kk] <= 0 or phiMVec[kk] >= 1:
            phiMVec[kk] = np.random.normal(0.192, 0.189)
            # phiMVec[kk] = np.random.rand()
        while dRatio[kk] <= 0 or dRatio[kk] >= 1:
            dRatio[kk] = np.random.normal(0.128, 0.021)
    
    # Experiment 20-25
    for kk in range(n//4):
        rWVec[kk+n//4] = np.random.normal(132/24, 25.5/576)
        # rWVec[kk+n//4] = np.random.rand()*4 + 4
        rMVec[kk+n//4] = np.random.normal(180/24, 55/576)
        while phiWVec[kk+n//4] <= 0 or phiWVec[kk+n//4] >= 1:
            phiWVec[kk+n//4] = np.random.normal(0.818, 0.255)
            # phiWVec[kk+n//4] = np.random.rand()
        while phiMVec[kk+n//4] <= 0 or phiMVec[kk+n//4] >= 1:
            phiMVec[kk+n//4] = np.random.normal(0.644, 0.55)
            # phiMVec[kk+n//4] = np.random.rand()
        while dRatio[kk+n//4] <= 0 or dRatio[kk+n//4] >= 1:
            dRatio[kk+n//4] = np.random.normal(0.308, 0.0404)

    # Experiment 25-20
    for kk in range(n//4):
        rWVec[kk+n//2] = np.random.normal(49/24, 23/576)
        # rWVec[kk+n//2] = np.random.rand()*4 + 0.5
        rMVec[kk+n//2] = np.random.normal(21.5/24, 37.3/576)
        while phiWVec[kk+n//2] <= 0 or phiWVec[kk+n//2] >= 1:
            phiWVec[kk+n//2] = np.random.normal(0.32, 0.23)
            # phiWVec[kk+n//2] = np.random.rand()
        while phiMVec[kk+n//2] <= 0 or phiMVec[kk+n//2] >= 1:
            phiMVec[kk+n//2] = np.random.normal(0.521, 0.0373)
            # phiMVec[kk+n//2] = np.random.rand()
        while dRatio[kk+n//2] <= 0 or dRatio[kk+n//2] >= 1:
            dRatio[kk+n//2] = np.random.normal(0.0254, 0.0121)
    
    # Experiment 25-25
    for kk in range(n//4):
        rWVec[kk+3*n//4] = np.random.normal(141/24, 41.8/576)
        # rWVec[kk+3*n//4] = np.random.rand()*4 + 4
        rMVec[kk+3*n//4] = np.random.normal(79/24, 26/576)
        while phiWVec[kk+3*n//4] <= 0 or phiWVec[kk+3*n//4] >= 1:
            phiWVec[kk+3*n//4] = np.random.normal(0.896, 0.418)
            # phiWVec[kk+3*n//4] = np.random.rand()
        while phiMVec[kk+3*n//4] <= 0 or phiMVec[kk+3*n//4] >= 1:
            phiMVec[kk+3*n//4] = np.random.normal(0.403, 0.26)
            # phiMVec[kk+3*n//4] = np.random.rand()
        while dRatio[kk+3*n//4] <= 0 or dRatio[kk+3*n//4] >= 1:
            dRatio[kk+3*n//4] = np.random.normal(0.198, 0.0591)

    # Convert these to the ND versions
    alphaWVec = -1/24*np.log(dRatio)
    gammaVec = np.random.rand(4*n//4)*1e-2

    # Loop through the n variable
    for kk in tqdm(range(4*(n//4)), leave = False, ncols = 50):
        
        # Create an empty dictionary for this repeat and populate with the
        # random parameters
        dfDict = {}
        if kk < n//4:
            dfDict['Experiment'] = ['20-20']
        elif kk < n//2:
            dfDict['Experiment'] = ['20-25']
        elif kk < (3*n)//4:
            dfDict['Experiment'] = ['25-20']
        else:
            dfDict['Experiment'] = ['25-25']
        dfDict['rW'] = [rWVec[kk]]
        dfDict['rM'] = [rMVec[kk]]
        dfDict['phiW'] = [phiWVec[kk]]
        dfDict['phiM'] = [phiMVec[kk]]
        dfDict['sigmaW'] = [alphaWVec[kk]]
        dfDict['gamma'] = [gammaVec[kk]]

        # Update pars
        pars['growthW'] = rWVec[kk]
        pars['growthM'] = rMVec[kk]
        pars['fecunW'] = phiWVec[kk]
        pars['fecunM'] = phiMVec[kk]
        pars['virW'] = alphaWVec[kk]
        pars['pathDil'] = gammaVec[kk]

        # Run a simulation
        sim = dimensionalSim(pars)
        sim.solveSystem()

        # Add to the dictionary with the outputs
        dfDict['SW'] = [sim.SS[0] if np.abs(sim.SS[0])>1e-3 else 0.0] 
        dfDict['IW'] = [sim.SS[1] if np.abs(sim.SS[1])>1e-3 else 0.0] 
        dfDict['SM'] = [sim.SS[2] if np.abs(sim.SS[2])>1e-3 else 0.0] 
        dfDict['IM'] = [sim.SS[3] if np.abs(sim.SS[3])>1e-3 else 0.0] 
        dfDict['N'] = [np.sum(sim.SS[:4]) if np.abs(np.sum(sim.SS[:4]))>1e-3 \
                    else 0.0]
        dfDict['P'] = [sim.SS[4] if np.abs(sim.SS[4])>1e-3 else 0.0] 

        # Append this repeat to the wider dataFrame
        df1 = pd.DataFrame().from_dict(dfDict)
        df = pd.concat([df, df1], ignore_index=True)

    # Save if necessary
    if save:

        # Generate a save name
        saveName = LU.createSaveName(f'DEXP')
        
        # Create the two directories in ./Data and ./Figures
        os.mkdir(f'./Data/{saveName}/')
        
        # Create a parameters file and store it in the Data directory
        del pars['growthW']
        del pars['growthM']
        del pars['fecunW']
        del pars['fecunM']
        del pars['virW']
        del pars['pathDil']

        with open(f'./Data/{saveName}/' + 'parameters.txt', 'w') as file:
            file.write('Parameters:\n')
            file.write('----------\n')
            file.write('\n')
            for key, val in pars.items():
                file.write(f'{key: <13}: {val:.3f}\n')
        
        # Save the dataframe
        df.to_csv(f'./Data/{saveName}/dataset.csv', index=False)

        # Add an entry to the lookup table
        LU.addEntrytoLookup(saveName, pars, 'DEXP')

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

def plotData(folderName, experiment, save = False):
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
    dfAll = pd.read_csv(f'./Data/{folderName}/dataset.csv')

    # Generate some new columns in the dataframe from old columns
    dfAll['XW'] = (dfAll['SW'] + dfAll['IW']) # Density of wildtype
    dfAll['XM'] = (dfAll['SM'] + dfAll['IM']) # Density of mutant
    dfAll['chi'] = (dfAll['XW'] - dfAll['XM'])/(dfAll['XW'] + dfAll['XM']) # Ratio of wildtype to mutant

    # Split the dataset by experiment
    if experiment == 'All':
        df = dfAll
    else:
        df = dfAll[dfAll['Experiment'] == experiment]

    # Generate a figure
    figWidth = 20
    figHeight = 20
    fsize = 20
    fig = plt.figure(figsize=(figWidth, figHeight))
    axr = fig.add_subplot(221)
    axDil = fig.add_subplot(222)
    axChi = fig.add_subplot(223)
    axComp = fig.add_subplot(224)
    fig.suptitle(f"Experiment: {experiment}", fontsize = fsize)

    # Plot 1: Growth ratio vs chi
    axr.scatter(df['rM']/df['rW'], df['chi'], c=df['chi'], 
                    s=20, cmap='coolwarm_r', vmin=-1.0, vmax=1.0)
    axr.set_xlim([min(df['rM']/df['rW']), max(df['rM']/df['rW'])])
    axr.set_xlabel(r'$\lambda$', fontsize=fsize)
    axr.set_xticks([min(df['rM']/df['rW']), max(df['rM']/df['rW'])])
    axr.set_xticklabels([f"{min(df['rM']/df['rW']):.2f}", 
                         f"{max(df['rM']/df['rW']):.2f}"], fontsize=fsize)
    axr.set_ylim([-1, 1])
    axr.set_ylabel(r'Susceptible genotype ratio, $\chi$', fontsize=fsize)
    axr.set_yticks([-1, 0, 1])
    axr.set_yticklabels(['-1', '0', '1'], fontsize=fsize)
    axr.text(min(df['rM']/df['rW']), 1.05, 'A.', fontsize=fsize)

    # Plot 2: Fecundity ratio and growth ratio vs chi
    axDil.scatter(df['gamma'], df['chi'], c=df['chi'], s=20, cmap='coolwarm_r',
                    vmin=-1.0, vmax=1.0)
    
    maxGamma = max(dfAll['gamma'])
    axDil.set_xlim([0, maxGamma])
    axDil.annotate('', xy=(0.95, -0.025), xycoords='axes fraction', 
                   xytext=(0.05, -0.025), 
                   arrowprops=dict(color='k'))
    axDil.text(maxGamma/2, -1.125, 
               "Increasing dilution effects", fontsize=fsize,
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

    # Create a secondary plot which looks at the losses in fecundity
    fig2 = plt.figure(figsize=(10, 5))
    ax21 = fig2.add_subplot(121)
    ax22 = fig2.add_subplot(122)

    x1, y1 = np.mgrid[0:1.01:0.01, 0:1.01:0.01]
    pos1 = np.dstack((x1, y1))
    x2, y2 = np.mgrid[0.5:8.01:0.01, 0:1.01:0.01]
    pos2 = np.dstack((x2, y2))
    if experiment == '20-20':
        mn1 = np.array([0.728, 0.192])
        cv1 = np.array([[0.233**2, 0],[0, 0.189**2]])
        mn2 = np.array([52.3/24, 0.728])
        cv2 = np.array([[(23.3/576)**2, 0],[0, 0.233**2]])
    if experiment == '20-25':
        mn1 = np.array([0.818, 0.633])
        cv1 = np.array([[0.255**2,0],[0,0.55**2]])
        mn2 = np.array([132/24, 0.818])
        cv2 = np.array([[(25.5/576)**2, 0],[0, 0.255**2]])
    if experiment == '25-20':
        mn1 = np.array([0.32, 0.521])
        cv1 = np.array([[0.23**2,0],[0,0.0373**2]])
        mn2 = np.array([49/24, 0.32])
        cv2 = np.array([[(23/576)**2, 0],[0, 0.23**2]])
    if experiment == '25-25':
        mn1 = np.array([0.896, 0.403])
        cv1 = np.array([[0.418**2,0],[0,0.26**2]])
        mn2 = np.array([141/24, 0.896])
        cv2 = np.array([[(41.8/576)**2, 0],[0, 0.418**2]])
            
    ax21.scatter(df['phiW'], df['phiM'], c=df['chi'], 
                s=20, cmap='coolwarm_r', vmin=-1.0, vmax=1.0)
    ax21.contour(x1, y1, sp.multivariate_normal(mn1, cv1).pdf(pos1), colors='k')
    ax21.set_xlabel(r'$\varphi_W$', fontsize=fsize)
    ax21.set_xticks([0, 0.5, 1])
    ax21.set_xticklabels(['0', '0.5', '1'], fontsize=fsize)
    ax21.set_ylabel(r'$\varphi_M$', fontsize=fsize)
    ax21.set_yticks([0, 0.5, 1])
    ax21.set_yticklabels(['0', '0.5', '1'], fontsize=fsize)

    ax22.scatter(df['rW'], df['phiW'], c=df['chi'], 
                s=20, cmap='coolwarm_r', vmin=-1.0, vmax=1.0)
    ax22.contour(x2, y2, sp.multivariate_normal(mn2, cv2).pdf(pos2),
                    colors='k')
    ax22.set_xlim((min(df['rW']), max(df['rW'])))
    ax22.set_xlabel(r'$r_W$', fontsize=fsize)
    ax22.set_xticks([min(df['rW']), max(df['rW'])])
    ax22.set_xticklabels([f"{min(df['rW']):.2f}", 
                          f"{max(df['rW']):.2f}"], fontsize=fsize)
    ax22.set_ylabel(r'$\varphi_W$', fontsize=fsize)
    ax22.set_yticks([0, 0.5, 1])
    ax22.set_yticklabels(['0', '0.5', '1'], fontsize=fsize)

    # Save if necessary
    if save:
        if not os.path.isdir(f'./Figures/{folderName}'):
            os.mkdir(f'./Figures/{folderName}')

        fig.savefig(f'./Figures/{folderName}/results-{experiment}.png', 
                bbox_inches='tight')
        fig2.savefig(f'./Figures/{folderName}/fecun-loss-{experiment}.png', 
                bbox_inches='tight')
        
    else:

        plt.show()

if __name__ == '__main__':

    # Testing
    if len(sys.argv) < 2:
        pars = initDim()
        df = generateData(pars, n=100)
        for ii in range(100):
            print(list(df.iloc[ii]))
    
    # Saving data
    elif sys.argv[1] == 'save':
        pars = initDim(tEco=50)
        generateData(pars, n = 2000, save = True)

    # Plotting in window
    elif sys.argv[1] == 'plot':
        plotData(sys.argv[2], '20-20')

    # Save plotting
    elif sys.argv[1] == 'saveplot':
        plotData(sys.argv[2], sys.argv[3], save=True)

    elif sys.argv[1] == 'saveplotall':
        for experiment in ['20-20', '20-25', '25-20', '25-25']:
            plotData(sys.argv[2], experiment, save=True)
