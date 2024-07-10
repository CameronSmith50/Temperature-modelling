# Functions which will be associated with a lookup table. Will be able to
# generate the file structure for the code, generate the lookup table and be
# able to add and remove entries from it.
#
# Author: Cameron Smith
# Date created: 28/11/2023
# Last modified: 09/04/2024

#%% Import any modules that are required

import os
import pandas as pd
from faker import Faker
import json
import numpy as np
import sys
import datetime
from copy import deepcopy
import shutil

#%% Functions

# Creation of a lookup table
def makeLookup(colNames, recreate=False):
    '''
    Code which will create a lookup table with column names given in 
    colnames.
    If the file exists, it can be recreated if recreates is set to True.
    '''

    if not os.path.isfile('./lookuptable.csv') or recreate == True:
        df = pd.DataFrame(columns=colNames)
        df.to_csv('./lookupTable.csv', index=False)

def addEntrytoLookup(saveName, pars, simType, data=True, figure=True):
    '''
    Adds entries to the lookup table. Requires the pars, which is a
    dictionary with the colNames from makeLookup as the keys, the save name for
    the files and the tyoe of simulation (string). Values for pars
    should be in a list. If data/figure are True, we have an entry in the
    Data/Figures directory.
    '''

    # Check if the lookup table exists
    if not os.path.isfile('./lookupTable.csv'):
        print('Lookup table does not yet exist. Please run "Lookup.py" in '
              'order to initialise it.')
        return
    
    # Create a copy of pars
    dataEntry = deepcopy(pars)
    
    # Add in the filename, time-stamp and which folders have the data
    dataEntry['Folder'] = saveName
    dataEntry['Data'] = data
    dataEntry['Figure'] = figure
    dataEntry['SimType'] = simType
    dt = datetime.datetime.now()
    dataEntry['Date'] = f'{dt.year}-{dt.month}-{dt.day}'
    dataEntry['Time'] = f'{dt.hour}:{dt.minute}'

    # Check if the columnames match with the columns of the dataframe
    df = pd.read_csv('./lookupTable.csv', index_col=False)
    dfCols = df.columns
    dataCols = list(dataEntry.keys())
    if not all(x in dfCols for x in dataCols):
        print("The inputted data columns don't match the lookup table")
        return

    # Loop through dataEntry and turn each item into a string
    for key, val in dataEntry.items():
        dataEntry[key] = [val]

    # For anything missing, add a dash
    for key in dfCols:
        if key not in dataCols:
            dataEntry[key] = ['-']
    
    # Now add an entry to the dataframe
    df1 = pd.DataFrame(dataEntry)
    df = pd.concat([df, df1], ignore_index=True)
    df.to_csv('./lookupTable.csv', index=False)

def removeEntryFromLookup(fileToRemove):
    '''
    Code to remove a data directory and its associated entry in the lookup
    table and any results associated with it. Will ask the user if they are
    sure they want to delete it before doing so. The fileToRemove is a string
    which is the three word identifier for the directories
    '''
    
    # Firstly, check if the file directory exists
    dataFileToRemove = f'./Data/{fileToRemove}/'
    if not os.path.isdir(dataFileToRemove):
        print(f'The directory "{dataFileToRemove}" does not exist, and so cannot '
              'be deleted. Have you remembered to end your directory with a '
              '"/"?')
        return
        
    # Otherwise, find the corresponding element in the lookup table
    print('')
    df = pd.read_csv('./lookupTable.csv')
    indToRem = df[df['Folder'] == fileToRemove].index[0]

    # Ask the user if they are sure they want to delete this file
    delete = input('Are you sure you want to delete the directories '
                   f'associated with "{fileToRemove}"? This will be a ' 'permanent change that cannot be undone. [y/n]: ')
    if delete == 'y':
        print('')
        print(f'Removing data directory "{dataFileToRemove}"...')
        shutil.rmtree(dataFileToRemove)
        print('DONE')
        print('Checking for associated results and deleting...')
        resultPath = dataFileToRemove.split('/')
        resultPath[1] = 'Figures'
        resultPath = '/'.join(resultPath)
        if os.path.isdir(resultPath):
            shutil.rmtree(resultPath)
            print('Results found and deleted')
        else:
            print('No figures directory found')
        print('Removing entry from the lookup table...')
        df = df.drop(indToRem)
        df.to_csv('./lookupTable.csv', index=False)
        print('DONE')
    else:
        print('')
        print('Files not deleted')

def createSaveName(prefix=None):
    '''
    This will use the words.json file to generate a unique name for a file
    directory. If prefix is a string, this will be prefixed to the name.
    '''

    # Check if the words.json file exists
    if not os.path.isfile('./words.json'):
        print('words.json does not yet exist. Please run "Lookup.py" in '
              'order to initialise it.')
        return
    
    # If it does, we choose three words at random with replacement
    with open("words.json", 'r') as file:
        wordList = json.load(file)

    # Reset prefix if none
    if prefix == None:
        prefix = ''
    else:
        prefix = prefix + '-'

    # Coose three integers at random between 0 and len(wordList)
    nmax = len(wordList)
    wordInd = [np.random.randint(0, nmax) for _ in range(3)]
    for ii, ind in enumerate(wordInd):
        wordInd[ii] = wordList[ind]

    return(f"{prefix}{'-'.join(wordInd)}")

#%% Initialisation
if __name__ == '__main__':
    
    if len(sys.argv) < 2:

        # Create the file directories if they don't exist
        if not os.path.isdir('./Data/'):
            os.makedirs('./Data/')

        if not os.path.isdir('./Figures/'):
            os.makedirs('./Figures/')

        # Create the lookup table
        colnames = [
            'Folder',
            'Data',
            'Figure',
            'SimType',
            'Date',
            'Time',
            'tEco',
            'growthW',
            'growthM',
            'fecunW',
            'fecunM',
            'transmission',
            'natMort',
            'virW',
            'virM',
            'comp',
            'pathShed',
            'pathDecay',
            'pathLoss',
            'pathDil',
            'pathGrow',
            'pathCap',
            'growthRatio',
            'natMortND',
            'virNDW',
            'virNDM',
            'pathShedND',
            'pathDecayND',
            'pathLossND',
            'pathDilND',
            'pathGrowND',
            'pathCapND'
            ]
        makeLookup(colNames=colnames)

        # Generate a list of 5000 words
        if not os.path.isfile('./words.json'):
            fake = Faker()
            words = [fake.word() for _ in range(5000)]
            with open('./words.json', 'w') as file:
                json.dump(list(set(words)), file)

        # Create a batch script for deletion. This will work by typing delFile
        # [fileToRemove]. the file to remove should
        # be the unique identifier only, not the entire directory (e.g.
        # ECO-single-strength-person)
        with open('delFile.bat', 'w') as file:
            file.write('python Lookup.py removeEntryFromLookup %1')

    # Code to be able to delete a directory from the command line. Should so
    # this using the delFile batch script above
    if len(sys.argv) > 2:
        globals()[sys.argv[1]](sys.argv[2])