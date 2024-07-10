# README

## Basic information

This is the accompanying code for the modelling section of the manuscript
"Warming during different life stages has distinct impacts on host resistance evolution" by Jingdi Li, Cameron Smith, Jinlin Chen and Kayla King. It is maintained by Cameron Smith
(cameronsmith50@outlook.com)

## Requirements and file structure

The code was written using Python 3.9.7. The following packages are required to
execute the code:
- `numpy`
- `matplotlib`
- `os`
- `sys`
- `json`
- `pandas`
- `faker`
- `copy`
- `shutil`
- `datetime`
- `tqdm`

There is no explicit file structure required. This will be generated when
needed.

## Code execution

To run the code, we firstly need to run the `Lookup.py` code by typing the
following while the code folder is the current working directory:

```sh
python Lookup.py
```
This will create the following file structure within the working directory
```sh
your_project_name/ <-- cd
|-Data/
|-Results/
.gitignore
delFile.bat
ecoSimDim.py
ecoSimND.py
experiments_v1.py
experiments_v2.py
Initialisation.py
Lookup.py
lookupTable.csv
paramSweep.py
README.md
words.json
```

The file `lookupTable.csv` will store all of the parameters for any saved
dataset, together with a time-stamp for when it was generated. The columns are
defined in `Lookup.py` and should match all of the parameters stored in the
`Initialisation.py` functions.

### Initialisation.py

`Initialisation.py` is a file contiaining all of teh default parameters. There
are two functions in this file, one for the dimensional system and one for the
non-dimensional system. Any paramters added to this file must be added to the
list in `Lookup.py` for saving to continue. Once added, Lookup needs to be rerun
as above.

If the user wishes to check that all of the parameters are correct, they may
execute the code 
```sh
python Lookup.py > Default_parameters.txt
```

which will create a .txt file with all of the default parameters in the root
directory. Omit anything including and after the ">" to print this to the
console instead.

### ecoSim*.py

The two functions `ecoSimDim.py` and `ecoSimND.py` act in the same way, but use
the dimensional and non-dimensional parameters respectively. To change the RHS
function of the ODEs, you can change the `RHSFunction` function, but if new
parameters are added then they need to be defined in `__init__`,
`Initialisation.py` and `Lookup.py`.

To run the code without saving any output, we can run 

```sh
python ecoSim*.py
```

where the * is replaced by either "Dim" or "ND". If the output is to be saved,
we execute the following:

```sh
python ecoSim*.py save
```

### paramSweep.py

`paramSweep` is the code used to generate the randomly parameterised plots seen
in the paper. To run this code we have four options:

#### Generate a small dataset and print to .txt file

To generate a small dataset and print to a txt file, we run the following code
in the command window

```sh
python paramSweep.py > test.txt
```

This should only be done as a test to ensure the code is running. By default
this runs 100 repeats.

#### Save a dataset

To save a much larger dataset which will be used to be processed later on, we
use the following code

```sh
python paramSweep.py save
```

This will save a dataset in a randomly named folder in the ./Data/ direcory of
the form `SCAT-word1-word2-word3` This folder will contain:
- `dataset.csv`: This is a file which has all the information the code requires
  to plot, including the randomly initialised parameters, steady states etc.
- `parameters.txt`: Contains all of the fixed parameters.

#### Plot a saved dataset without saving

In order to plot a dataset that has been saved as above, we type

```sh
python paramSweep plot SCAT-word1-word2-word3
```

Again, this should be done when changing any figures so that you don't have to
keep overwriting saved figures.

### Plot and save using a saved dataset

To generate a plot for the dataset which is saved, we execute the following code

```sh
python paramSweep.py saveplot SCAT-word1-word2-word3
```

where a new folder is created `./Figures/SCAT-word1-word2-word3/` with the
figure inside.

### experiments_v*.py

These two python scripts act in a similar way to the `paramSweep` file, but use
experimental distributions instead of completely random parameter sets.
`experiments_v1.py` uses the non-dimensional parameter sets while `experiments_v2.py`
uses the dimensional set. The code is executed in a similar way to
`paramSweep.py`, but with the prefixes "EXP" or "DEXP" respectively instead of
"SCAT". If using plot, the 20-20 experiment will plot by default. If using
saveplot, an experiment must be chosen to be plotted, e.g. 20-20, 20-25, 25-20
or 25-25 as a final argument. For example the following would plot the results
from the 20-25 experiment stored in DEXP-word1-word2-word3:

```sh
python experiment_v2.py saveplot DEXP-word1-word2-word3 20-25
```

Finally, a `saveplotall` keyword can be used instead, such as below:

```sh
python experiment_v2.py saveplotall DEXP-word1-word2-word3
```

which will save all experiments.

## Other files

`delFile.bat` is used to delete any saved directories. It is used in the
following way:

```sh
delFile SCAT-word1-word2-word3
```
which would delete any data or figure directories with that name, as well as its
entry from the file `lookupTable.csv`.