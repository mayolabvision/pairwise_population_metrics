function [ loading_similarity, L_orth ] = compute_load_sim( fa_model )
%
% computes the loading similarity of each dimension
%
% Input:
%   fa_model - (struct), factor analysis model parameters
%
% Output:
%   loading_similarity - (list), loading similarities of each dimension's
%                        loading vector
%   L_orth              - (n_neurons x z_dim), orthonormalized per-neuron
%                        loadings (same matrix loading_similarity is
%                        computed from) -- e.g. abs(L_orth(:,1)) is each
%                        neuron's loading magnitude onto the dominant
%                        shared dimension
%
% @ Akash Umakantha, 2021. See https://www.biorxiv.org/content/10.1101/2020.12.04.383604v1
%

    % return empty if 'L' does not exist
    if ~isfield(fa_model,'L')
        loading_similarity=[];
        L_orth=[];
        return
    end

    % extract some basic info/params
    L = fa_model.L;
    [n_neurons,z_dim] = size(L);
    
    % orthonormalize L
    shared_cov = L*L';
    [L_orth,~,~] = svd(shared_cov);
    L_orth = L_orth(:,1:z_dim);
    
    % compute loading similarity of each dimension
    loading_similarity = (1/n_neurons - var(L_orth,1,1)).*n_neurons;
    
end
