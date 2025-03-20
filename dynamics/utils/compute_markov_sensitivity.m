function [sg_rel, h_rel] = compute_markov_sensitivity(agg_tm, agg_spectral_gap, agg_h_rate, num_samples, lowerBound, upperBound)
% compute_markov_sensitivity - Computes the sensitivity of Markov dynamics to perturbations.
%
% Syntax:
%   [sg_rel, h_rel] = compute_markov_sensitivity(agg_tm, agg_spectral_gap, agg_h_rate, num_samples, lowerBound, upperBound)
%
% Inputs:
%   agg_tm           - Aggregate transition matrix (n x n) of the Markov chain.
%   agg_spectral_gap - Baseline spectral gap computed from agg_tm.
%   agg_h_rate       - Baseline entropy rate (in bits) computed from agg_tm.
%   num_samples      - Number of surrogate samples to generate (perturbations).
%   lowerBound       - Lower bound for the perturbation magnitude (epsilon).
%   upperBound       - Upper bound for the perturbation magnitude (epsilon).
%
% Outputs:
%   sg_rel - A 3D array (num_samples x n x n) containing the relative change in the 
%            spectral gap (in percentage) for each perturbation.
%   h_rel  - A 3D array (num_samples x n x n) containing the relative change in the 
%            entropy rate (in percentage) for each perturbation.
%
% Description:
%   For each surrogate sample, a random epsilon is generated (in log space between lowerBound
%   and upperBound) and a random initial distribution is created. For every entry (i,j) in the
%   aggregate transition matrix, epsilon is added to that entry, and then the i-th row is renormalized.
%   The function computes the perturbed stationary distribution, spectral gap, and entropy rate using
%   the helper functions compute_markov_dynamics and markovEntropyRate. The relative changes with respect
%   to the original (baseline) values are stored in sg_rel and h_rel.
%
% Example:
%   [sg_rel, h_rel] = compute_markov_sensitivity(agg_tm, agg_spectral_gap, agg_h_rate, 10, 1e-4, 1e-2);
%
% Author: [Sir-Lord]
% Date: [Today's Date]

    % Determine the number of states (clusters)
    num_clusters = size(agg_tm, 1);
    
    % Preallocate arrays to store relative changes for spectral gap and entropy rate.
    sg_rel = zeros(num_samples, num_clusters, num_clusters);
    h_rel = zeros(num_samples, num_clusters, num_clusters);
    
    % Compute logarithms of the lower and upper bounds for epsilon.
    logMin = log(lowerBound);
    logMax = log(upperBound);
    
    % Loop over the number of surrogate samples.
    for s = 1:num_samples
        fprintf("Sample: (%d/%d)\n", s, num_samples)
        % Generate a random epsilon in log space.
        eps = exp(logMin + (logMax - logMin) * rand);
        
        % Generate a random initial distribution (normalized).
        rnd_init = rand(num_clusters, 1);
        rnd_init = rnd_init / sum(rnd_init);
        
        % Loop over each entry (i,j) of the transition matrix.
        for idx_1 = 1:num_clusters
            for idx_2 = 1:num_clusters
                % Create a perturbed transition matrix by copying the original.
                agg_tm_pert = agg_tm;
                % Add epsilon to the (idx_1, idx_2) entry.
                agg_tm_pert(idx_1, idx_2) = agg_tm(idx_1, idx_2) + eps;
                % Renormalize the perturbed row to ensure it sums to 1.
                agg_tm_pert(idx_1, :) = agg_tm_pert(idx_1, :) / sum(agg_tm_pert(idx_1, :));
                
                % Compute the perturbed Markov dynamics.
                [agg_stat_dist_pert, agg_spectral_gap_pert, ~, agg_tv_steps_pert] = ...
                    compute_markov_dynamics(agg_tm_pert, 'mu_initial', rnd_init);
                
                % Calculate relative change in spectral gap (in percentage).
                sg_rel(s, idx_1, idx_2) = 100 * ((agg_spectral_gap_pert - agg_spectral_gap) / agg_spectral_gap);
                
                % Compute the perturbed entropy rate.
                agg_h_rate_pert = markovEntropyRate(agg_stat_dist_pert, agg_tm_pert);
                % Calculate relative change in entropy rate (in percentage).
                h_rel(s, idx_1, idx_2) = 100 * ((agg_h_rate_pert - agg_h_rate) / agg_h_rate);
            end
        end
    end
end
