function [ls, psv, d, eigspec, lsN, psvN] = compute_population_metrics(X, zDims)
%
% fits factor analysis (using crossvalidaiton to select a dimensionality)
% and then computes population metrics and shared eigenspectrum
%
% Input:
%   X     - (matrix of size n_neurons x n_trials), spike count matrix
%   zDims - (list of integers), dimensionalities to try for FA crossvalidation
%
% Output:
%   ls      - (list), loading similarity of each fitted dimension
%   psv     - (scalar), percent shared variance, averaged across neurons
%   d       - (scalar), shared dimensionality (dshared)
%   eigspec - (list), eigenspectrum of shared covariance matrix
%   lsN     - (n_neurons x z_dim), orthonormalized per-neuron loadings
%             ls was computed from
%   psvN    - (n_neurons x 1), percent shared variance per neuron
%
% @ Akash Umakantha, 2021. See https://www.biorxiv.org/content/10.1101/2020.12.04.383604v1
%

    if nargin<2
        zDims = 0:10;
    end

    % compute cross-validated LL for each dimensionality
    dims = crossvalidate_fa(X,'zDimList',zDims);
    
    % find the model that gives highest cross-validated LL
    [~,max_idx] = max([dims.sumLL]);
    params = dims(max_idx).estParams;
    
    % compute loading similarity, %sv, dshared, and shared eigenspectrum
    [ls, lsN] = compute_load_sim(params);
    [psv, psvN] = compute_perc_shared(params);
    d = compute_dshared(params);
    eigspec = compute_shared_eigspec(params);

end

