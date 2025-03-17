# Code by Melinda Kleczynski 

# Finalized February 27, 2025 

# Input a file path, values of epsilon, values of sigma 
# Read a csv file with 3d spatial coordinates, chain labels, and residue numbering of atoms 
# Return a list of GCCD matrices, one for each value of sigma

import numpy as np 
import oat_python as oat
import pandas as pd 
from scipy.sparse import csr_matrix 
from scipy.spatial.distance import pdist, squareform
from scipy.special import ndtr  # cdf of normal distribution 

def smooth_betti(epsilon_scalar, births_vec, deaths_vec, sigma_scalar):

    scaled_births = (epsilon_scalar - births_vec)/sigma_scalar 
    scaled_deaths = (epsilon_scalar - deaths_vec)/sigma_scalar 
    
    return np.sum(ndtr(scaled_births) * (1 - ndtr(scaled_deaths)))

def gccd(atom_data_fpath, epsilons, sigmas):

    # set up data
    frame_data = pd.read_csv(atom_data_fpath) 
    n_atoms = len(frame_data)
    atom_chains = np.array(frame_data.chain)
    atom_residues = np.array(frame_data.residue)
    atom_coords = np.array(frame_data[['x', 'y', 'z']])

    # determine maximum residue distance
    max_inchain_resdist = 0
    chains = np.unique(atom_chains)
    for chain in chains:
        residues_in_chain = frame_data[frame_data.chain == chain].residue
        inchain_resdist = np.max(residues_in_chain) - np.min(residues_in_chain)
        max_inchain_resdist = np.maximum(inchain_resdist, max_inchain_resdist)
    if len(chains) > 1:
        max_resdist = max_inchain_resdist + 1  # last column for features requiring both chains 
    else:
        max_resdist = max_inchain_resdist 

    # compute pairwise residue distances 
    res_pdist = np.zeros((n_atoms, n_atoms), int)
    for i in range(n_atoms):
        for j in range(i):
            if atom_chains[i] != atom_chains[j]:
                res_pdist[i, j] = max_resdist 
            else:
                res_pdist[i, j] = np.abs(atom_residues[i] - atom_residues[j])
    res_pdist += np.transpose(res_pdist)

    # compute pairwise spatial distances 
    end_filtration = epsilons[-1] + 0.001 
    spatial_pdist = squareform(pdist(atom_coords, metric = 'euclidean'))
    spatial_pdist = spatial_pdist * (spatial_pdist <= end_filtration)

    # compute Gaussian CROCKER matrices
    gauss_crocks = [np.zeros((len(epsilons), max_resdist + 1)) for s in range(len(sigmas))] # will leave the first column as all zeros for when we take the difference 
    for resdist_threshold in range(1, max_resdist + 1): 

        # format pairwise distance matrix 
        filtration_pdist = spatial_pdist * (res_pdist <= resdist_threshold)
        formatted_pdist = csr_matrix(filtration_pdist) 
        formatted_pdist[[i for i in range(n_atoms)], [i for i in range(n_atoms)]] = [0 for i in range(n_atoms)]

        # run tda 
        boundary = oat.rust.FactoredBoundaryMatrixVr(dissimilarity_matrix = formatted_pdist, homology_dimension_max = 1)  
        ph = boundary.homology(return_cycle_representatives = False, return_bounding_chains = False) 

        # extract births and deaths 
        ph_dim1 = ph[ph.dimension == 1]
        births = np.array(ph_dim1.birth)
        deaths = np.array(ph_dim1.death) 
        deaths = np.nan_to_num(deaths, posinf = end_filtration)

        # compute smooth Betti curves 
        if len(deaths) > 0:
            for s in range(len(sigmas)):
                gauss_crocks[s][:, resdist_threshold] = [smooth_betti(epsilon, births, deaths, sigmas[s]) for epsilon in epsilons]  

    # compute GCCDs
    gccds = [gauss_crocks[s][:, 1:] - gauss_crocks[s][:, 0:-1] for s in range(len(sigmas))]

    return gccds 