%% Clear Environment and Set Paths
clc;clear;
addpath(genpath('./'));

%% Data Initialization
% Define parameters for simulated fMRI time courses.
num_subjects   = 50;         % Number of subjects
num_timepoints = 150;         % Number of timepoints per subject
num_components = 5;         % Number of components per subject
num_features   = nchoosek(num_components, 2); % Number of unique component pairs
Tr             = 0.72;       % Repetition time (s)
Fs             = 1/Tr;       % Sampling frequency (Hz)
nyquist        = Fs/2;       % Nyquist frequency
band           = [0.01 0.198];  % Frequency band for post-processing
cutoff         = [0.6 1.4];     % Filter cutoff parameters

% Generate example post-processed fMRI timecourses for each session.
% The post_processing function is assumed to be available.
subTcs = randn(num_subjects, num_timepoints, num_components);
subTcs = post_processing(subTcs, Tr, band, cutoff);

%% DTW Parameters
% Define a range of gamma values for the DTW computation.
gamma_vec = 1.5;
num_gammas = length(gamma_vec);

%% DTW and nDTW Computation across all subjects and gammas.
[~, ~, full_tr_DTW, full_d_DTW] = compute_subject_dtw_framework(subTcs, gamma_vec, Tr, band);
full_tr_DTW = squeeze(full_tr_DTW(1, :, :, :));
full_d_DTW = squeeze(full_d_DTW(1, :, :, :));
%% Groupings
% Define group labels: schizophrenia subjects (1) and controls (0).
sz_label = 1;
cn_label = 0;
group = randi([cn_label, sz_label], num_subjects, 1);  % Random binary group assignment

% Separate subjects into two groups based on cluster index
sz_idx = find(group == sz_label);
cn_idx = find(group == cn_label);

% Count the number of subjects in each group
num_sz = length(sz_idx);
num_cn = length(cn_idx);
%%  Perform K-means
% Concatenate subjects and time of the trFNC to get a 2D matrix
tr_DTW_k = reshape(full_tr_DTW, num_subjects*num_timepoints, num_features);

% Define the number of clusters
num_clusters = 3;
replications = 20;

% Perform k-means
% [cluster_idx, cluster_c, cluster_sumd, cluster_D] = calculate_kmeans(we_k, cluster_num, 'cityblock', replications);
% cluster_idx = reshape(cluster_idx, subjects_cnt, []); 

% Generate random cluster indices for trDTW to mimic kmeans index output
cluster_idx = randi([1 num_clusters], num_subjects * num_timepoints, 1);

% Reshape into subjects x timepoints
cluster_idx = reshape(cluster_idx, num_subjects, num_timepoints);

%% Extract groups from metrics
% Extract trFNC data for each group
tr_DTW_sz = full_tr_DTW(sz_idx, :, :);
tr_DTW_cn = full_tr_DTW(cn_idx, :, :);

% Reshape cluster indices for each group
cluster_idx_sz = cluster_idx(sz_idx, :);
cluster_idx_cn = cluster_idx(cn_idx, :);

%% Extract subject states for each group
% Calculate group centroids for each group
sub_centroids_sz = group_sub_centroids(tr_DTW_sz, cluster_idx_sz, num_clusters);
% sub_centroids_cn = group_sub_centroids(tr_DTW_cn, cluster_idx_cn, num_clusters);

%% Analyze groups' cluster dynamics
% Compute cluster analysis metrics for each subject in each group
[mdt, fr, tm, tc] = trFNC_cluster_analysis(cluster_idx, num_clusters);
[mdt_sz, fr_sz, tm_sz, tc_sz] = trFNC_cluster_analysis(cluster_idx_sz, num_clusters);
[mdt_cn, fr_cn, tm_cn, tc_cn] = trFNC_cluster_analysis(cluster_idx_cn, num_clusters);

% Perform statistucan tests on mean dwell time (mdt), fractional occupancy
% (fr), transition matrix (tm) and transition count (tc)

%% Markov chain analysis
%Aggregate all subjects transition counts into aggregated transition
%matrix
agg_tm = aggregate_transitions(tc);
tol = 1e-3; % tolerance level

%Check if the transition matrix is ergodic
isErgodic = checkIsErgodic(agg_tm);
if isErgodic
    [agg_stat_dist, agg_spectral_gap, agg_tv_distance, agg_tv_steps] = compute_markov_dynamics(agg_tm, 'mu_initial', [1 0 0]', 'tol', tol);
    agg_h_rate = markovEntropyRate(agg_stat_dist, agg_tm);
end

%% Subject-wise markov chain analysis
[stat_dist_cn, spec_gap_cn, tv_steps_cn, hr_cn] = compute_sub_markov_analysis(tm_cn, tol);
[stat_dist_sz, spec_gap_sz, tv_steps_sz, hr_sz] = compute_sub_markov_analysis(tm_sz, tol);


%% Markov chain perturbation analysis
num_samples = 10000;
lowerBound = 1e-3;
upperBound = 1e-1;
[sg_rel, h_rel] = compute_markov_sensitivity(agg_tm, agg_spectral_gap, agg_h_rate, num_samples, lowerBound, upperBound);

%% 
function [sub_stat_dist, sub_spec_gap, sub_tv_steps, sub_hr] = compute_sub_markov_analysis(tm, tol)
% compute_sub_markov_analysis - Computes markov analysis and entropy rates for a group of subjects.
%
% Syntax:
%   [sub_stat_dist, sub_spec_gap, sub_tv_steps, sub_hr] = compute_sub_markov_analysis(tm, tol)
%
% Inputs:
%   tm  - A 3D array of transition matrices with dimensions:
%         [num_subjects x num_clusters x num_clusters], where each slice (along the first dimension)
%         represents the transition probability matrix for one subject.
%
%   tol - A tolerance value (e.g., 1e-10) to be used in the compute_state_vector function.
%
% Outputs:
%   sub_stat_dist - A matrix (num_subjects x num_clusters) where each row is the stationary
%                   distribution (state vector) computed for a subject.
%
%   sub_spec_gap  - A column vector (num_subjects x 1) of spectral gap values for each subject.
%
%   sub_tv_steps  - A column vector (num_subjects x 1) of total variation steps for each subject.
%
%   sub_hr        - A column vector (num_subjects x 1) of Markov chain entropy rates (in bits)
%                   computed for each subject.
%
% Description:
%   For each subject, the function extracts the transition probability matrix, checks whether
%   the Markov chain is ergodic and if each row sums to 1 (within tolerance). If both conditions
%   are satisfied, the function computes the stationary distribution, spectral gap, and total 
%   variation steps via compute_state_vector, and then calculates the entropy rate using 
%   markovEntropyRate. If the conditions are not met, the outputs for that subject remain NaN.
%
% Example:
%   [stat_dist, spec_gap, tv_steps, hr] = compute_sub_markov_analysis(tm, 1e-10);
%
% Author: [Sir-Lord]

    % Determine the number of subjects and clusters.
    num_subjects = size(tm, 1);
    num_clusters = size(tm, 2);
    
    % Preallocate output arrays.
    sub_stat_dist = nan(num_subjects, num_clusters);
    sub_spec_gap  = nan(num_subjects, 1);
    sub_tv_steps  = nan(num_subjects, 1);
    sub_hr        = nan(num_subjects, 1);
    
    % Loop over each subject.
    for sub = 1:num_subjects
        % Extract the transition matrix for the current subject.
        tm_sub = squeeze(tm(sub, :, :));
        
        % Check if the subject's Markov chain is ergodic.
        isErgodic = checkIsErgodic(tm_sub);
        
        % Compute row sums for validation.
        row_sums = sum(tm_sub, 2);
        
        % Process only if the chain is ergodic and all rows sum to 1 (within tolerance).
        if isErgodic && ~(any(abs(row_sums - 1) > 1e-10))
            % Compute stationary distribution, spectral gap, and total variation steps.
            [sub_stat_dist(sub, :), sub_spec_gap(sub), ~, sub_tv_steps(sub)] = ...
                compute_markov_dynamics(tm_sub, 'tol', tol);

            % Compute the entropy rate of the Markov chain using the stationary distribution.
            sub_hr(sub) = markovEntropyRate(sub_stat_dist(sub, :)', tm_sub);
        end
    end
end

