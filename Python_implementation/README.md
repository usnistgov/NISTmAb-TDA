# Python implementation

The main function file is gccd.py. Any notebook containing "from gccd import..." should have gccd.py saved in the same folder as the notebook to run.

## Python_versions.md

Versions used for the Python implementation at the time of writing.

In particular, note that Open Applied Topology (OAT) was recently released at the time of writing. Changes to the code may occur for later versions.

## call_betti_012.ipynb

Read csv files with 3d spatial coordinates, save concatenated 0d, 1d, and 2d normalized Gaussian Betti curves. 

## call_gccd.ipynb

Notebook for reading formatted data and saving GCCD matrices. 

## determine_filtration_range.ipynb

Determine when we can end the filtration / which rows to include in the GCCD matrices.

## example_pdb.ipynb

GCCD example on a different protein type, from PDB format.

## gccd.py

Main function file.

## test_accuracy.ipynb 

Classification test accuracy for three possible types of topological summaries: 
* GCCD matrices
* Gaussian Betti curves
* Concatenated normalized 0D, 1D, 2D Gaussian Betti curves

## time_gccd_fn.ipynb 

Time the GCCD function for a random selection of frames.

## tune_hyperparams.ipynb

Tune classification hyperparameters.
