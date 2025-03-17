# Julia implementation

See the Python implementation for hyperparameter tuning and classification.

## Julia_versions.md

Versions used for the Julia implementation at the time of writing.

In particular, note that the most recent release of Eirene (the library used for TDA computations) was in 2021. We were not able to use all of the most recent versions due to dependency requirements.

## call_gccd_Julia.ipynb

Read formatted data (csv files with 3d spatial coordinates, chain labels, and residue numbering of atoms), save GCCD matrices to csv files.

## gccd_Julia_functions.jl

Read a csv file with 3d spatial coordinates, chain labels, and residue numbering of atoms. Return a GCCD matrix for each value of sigma. 

## time_gccd_fn_Julia.ipynb

Time GCCD function.
