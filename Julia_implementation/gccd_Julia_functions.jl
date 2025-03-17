# Code by Melinda Kleczynski 

# Finalized March 5, 2025

# Read a csv file with 3d spatial coordinates, chain labels, and residue numbering of atoms 
# Return a GCCD matrix for each value of sigma 

using CSV 
using DataFrames
using Distances 
using Distributions 
using Eirene 

function smooth_betti(epsilon_scalar, barcode_2darray, sigma_scalar)

    n_bars = size(barcode_2darray, 1)

    return sum([cdf(Normal(barcode_2darray[bar_iter, 1], sigma_scalar), epsilon_scalar) *      # x (birth) ≤ ε
                (1 - cdf(Normal(barcode_2darray[bar_iter, 2], sigma_scalar), epsilon_scalar))  # y (death) > ε
                for bar_iter in 1:n_bars]) 

end 

function make_gccd_Julia(atom_data_fpath, epsilons, sigmas)

    # set up data 
    frame_data = CSV.read(atom_data_fpath, DataFrame)
    n_atoms = nrow(frame_data)
    atom_chains = frame_data.chain
    atom_residues = frame_data.residue
    atom_coords = hcat(frame_data.x, frame_data.y, frame_data.z)

    # determine maximum residue distance
    max_inchain_resdist = 0
    chains = unique(atom_chains)
    for chain in chains 
        residues_in_chain = atom_residues[atom_chains .== chain]
        inchain_resdist = maximum(residues_in_chain) - minimum(residues_in_chain)
        max_inchain_resdist = max(inchain_resdist, max_inchain_resdist)
    end 
    if length(chains) > 1
        max_resdist = max_inchain_resdist + 1  # last column for features requiring both chains 
    else
        max_resdist = max_inchain_resdist 
    end 

    # compute pairwise residue distances 
    res_pdist = zeros(Int, (n_atoms, n_atoms))
    for i in 1:n_atoms
        for j in 1:i-1 
            if atom_chains[i] != atom_chains[j]
                res_pdist[i, j] = max_resdist 
            else 
                res_pdist[i, j] = abs(atom_residues[i] - atom_residues[j])
            end 
        end
    end 
    res_pdist += transpose(res_pdist)

    # compute pairwise spatial distances 
    end_filtration = epsilons[end] + 0.001 
    spatial_pdist = pairwise(Euclidean(), atom_coords, dims = 1)

    # compute Gaussian CROCKER matrices
    gauss_crocks = [zeros((length(epsilons), max_resdist + 1)) for _ in 1:length(sigmas)] # will leave the first column as all zeros for when we take the difference

    for resdist_threshold in 1:max_resdist

        # give forbidden connections pairwise distances past the end of the filtration  
        filtration_pdist = spatial_pdist + (end_filtration + 5)*(res_pdist .> resdist_threshold)

        # run tda 
        eirene_output = eirene(filtration_pdist, model = "vr", maxrad = end_filtration, maxdim = 1) 
        bc1 = barcode(eirene_output, dim = 1)
        bc1 = replace(bc1, Inf => end_filtration)

        # compute smooth Betti curves
        if size(bc1)[1] > 0 
            for s in 1:length(sigmas) 
                gauss_crocks[s][:, resdist_threshold+1] = [smooth_betti(epsilon, bc1, sigmas[s]) for epsilon in epsilons]  
            end 
        end 

    end 

    # compute GCCDs
    gccds = [gauss_crocks[s][:, 2:end] - gauss_crocks[s][:, 1:end-1] for s in 1:length(sigmas)]

    return gccds 

end 