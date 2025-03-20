%% Clear Environment and Set Paths
clc; clear;
addpath(genpath('./'));

%% Data Initialization
% Define parameters for simulated fMRI time courses.
num_subjects   = 50;                          % Number of subjects
num_timepoints = 150;                         % Number of timepoints per subject
num_components = 5;                           % Number of components per subject
num_features   = nchoosek(num_components, 2);  % Number of unique component pairs
Tr             = 0.72;                        % Repetition time (s)
Fs             = 1/Tr;                        % Sampling frequency (Hz)
nyquist        = Fs/2;                        % Nyquist frequency

% Define processing bands and cutoffs for two pipelines.
band_f1    = [0.01 0.198];    % Frequency band for post-processing (f1)
band_f2    = [0.01 0.15];     % Frequency band for post-processing (f2)
cutoff_f1  = [0.6 1.4];       % Filter cutoff parameters (f1)
cutoff_f2  = [0.6 1.5];       % Filter cutoff parameters (f2)

% Generate example post-processed fMRI time courses.
% The post_processing function is assumed to be available.
subTcs = randn(num_subjects, num_timepoints, num_components);
subTcs_f1 = post_processing(subTcs, Tr, band_f1, cutoff_f1);
subTcs_f2 = post_processing(subTcs, Tr, band_f2, cutoff_f2);

%% DTW Parameters
% Define a range of gamma values for the DTW computation.
gamma_vec = 1;
num_gammas = length(gamma_vec);

%% DTW and nDTW Computation across All Subjects and Gammas
[full_dtw_d_f1, full_dtw_d_norm_f1] = compute_subject_dtw_framework(subTcs_f1, gamma_vec, Tr, band_f1);
[full_dtw_d_f2, full_dtw_d_norm_f2] = compute_subject_dtw_framework(subTcs_f2, gamma_vec, Tr, band_f2);

% For illustration, pick one gamma value (first one).
dtw_d_norm_f1 = squeeze(full_dtw_d_norm_f1(1, :, :));
dtw_d_norm_f2 = squeeze(full_dtw_d_norm_f2(1, :, :));

% Normalize the DTW for better comparison
dtw_d_norm_f1 = -zscore(dtw_d_norm_f1);
dtw_d_norm_f2 = -zscore(dtw_d_norm_f2);

%% Setting up Covariates (Group Differences)
% Define group labels: schizophrenia subjects (1) and controls (0).
sz_label = 1;
cn_label = 0;
group = randi([cn_label, sz_label], num_subjects, 1);  % Random binary group assignment

% Separate normalized DTW distances based on group membership.
dtw_d_norm_sz_f1 = dtw_d_norm_f1(group == sz_label, :);
dtw_d_norm_cn_f1 = dtw_d_norm_f1(group == cn_label, :);

dtw_d_norm_sz_f2 = dtw_d_norm_f2(group == sz_label, :);
dtw_d_norm_cn_f2 = dtw_d_norm_f2(group == cn_label, :);

%% Statistical Comparison between f1 and f2 for Each Feature
% Preallocate arrays for t-test p-values and effect sizes.
p_val_cn = zeros(1, num_features);
p_val_sz = zeros(1, num_features);
d_val_cn = zeros(1, num_features);
d_val_sz = zeros(1, num_features);

% Loop over each unique component pair (feature)
for f = 1:num_features
    % For controls: compare normalized DTW distances from pipeline f1 vs. f2.
    [~, p_val_cn(f)] = ttest(dtw_d_norm_cn_f1(:, f), dtw_d_norm_cn_f2(:, f));
    d_val_cn(f) = compute_effect_size(dtw_d_norm_cn_f1(:, f), dtw_d_norm_cn_f2(:, f)); 

    % For schizophrenia subjects: compare normalized DTW distances from pipeline f1 vs. f2.
    [~, p_val_sz(f)] = ttest(dtw_d_norm_sz_f1(:, f), dtw_d_norm_sz_f2(:, f));
    d_val_sz(f) = compute_effect_size(dtw_d_norm_sz_f1(:, f), dtw_d_norm_sz_f2(:, f)); 
end

% Apply FDR correction (Benjamini-Hochberg) to the p-values.
q_val_cn = mafdr(p_val_cn, "BHFDR", true);
q_val_sz = mafdr(p_val_sz, "BHFDR", true);

%% Function: Compute Effect Size (Cohen's d)
function d = compute_effect_size(x, y)
    % compute_effect_size - Computes Cohen's d effect size between two groups.
    %
    % Syntax:
    %   d = compute_effect_size(x, y)
    %
    % Inputs:
    %   x - A vector of measurements from group 1.
    %   y - A vector of measurements from group 2.
    %
    % Outputs:
    %   d - Cohen's d effect size.
    %
    % Note:
    %   The group with the lower label is taken as group1 and the group with the higher
    %   label is taken as group2.
    
    % Compute group means.
    mean1 = mean(x);
    mean2 = mean(y);

    % Compute group standard deviations.
    std1 = std(x);
    std2 = std(y);

    % Get sample sizes.
    n1 = length(x);
    n2 = length(y);

    % Compute the pooled standard deviation.
    pooled_std = sqrt(((n1 - 1)*std1^2 + (n2 - 1)*std2^2) / (n1 + n2 - 2));

    % Compute Cohen's d.
    d = (mean1 - mean2) / pooled_std;
end
