%% Clear Environment and Set Paths
clc;clear;
addpath(genpath('./'));

%% Data Initialization
% Define parameters for simulated fMRI time courses.
num_subjects  = 50;         % Number of subjects
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
gamma_vec = 0.25:0.25:3;
num_gammas = length(gamma_vec);

%% DTW and nDTW Computation across all subjects and gammas.
[full_dtw_d, full_dtw_d_norm] = compute_subject_dtw_framework(subTcs, gamma_vec, Tr, band);

%% Setting up group difference covariates
group = randi([0, 1], num_subjects, 1);      % Random binary group assignment (1: schizophrenia, 0: controls)

% Creating a table of covariates
age = randi([20, 80], num_subjects, 1);                % Random ages between 20 and 80
sex = randi([0, 1], num_subjects, 1);                  % Random binary sex assignment
site = randi([1, 3], num_subjects, 1);                 % Random site assignment with 3 possible sites
mfd = randn(num_subjects, 1);                          % Random mean frame displacement assignment
covariates = table(age, sex, site, mfd);    % Table of covariates

%% Compute group difference for each gamma value
% Preallocate output matrices.
dtw_d_q_val      = zeros(num_gammas, num_features); %q_val: FDR corrected p-value
dtw_d_p_val      = zeros(num_gammas, num_features); %p_val: Uncorrected p-value
dtw_d_t_val      = zeros(num_gammas, num_features); %t_val: t-statistic
dtw_d_d_val      = zeros(num_gammas, num_features); %d_val: cohen's d

dtw_d_norm_q_val = zeros(num_gammas, num_features); %q_val: FDR corrected p-value
dtw_d_norm_p_val = zeros(num_gammas, num_features); %p_val: Uncorrected p-value
dtw_d_norm_t_val = zeros(num_gammas, num_features); %t_val: t-statistic
dtw_d_norm_d_val = zeros(num_gammas, num_features); %d_val: cohen's d

% Loop over each gamma value.
for gamma_idx = 1:num_gammas
    % Display progress information in the command window.
    fprintf('Processing gamma: (%d/%d)\n', ...
        gamma_idx, num_gammas);
    % Extract the DTW features for the current gamma.
    dtw_d_vec      = squeeze(full_dtw_d(gamma_idx, :, :));      % [num_subjects x num_features]
    dtw_d_norm_vec = squeeze(full_dtw_d_norm(gamma_idx, :, :));   % [num_subjects x num_features]
    
    % Compute GLM statistics for vanilla DTW features.
    [q_values, p_values, t_values] = compute_glm_group_difference(dtw_d_vec, group, covariates);
    [d_values] = compute_cohens_d_effect(dtw_d_vec, group);
    dtw_d_q_val(gamma_idx, :) = q_values;
    dtw_d_p_val(gamma_idx, :) = p_values;
    dtw_d_t_val(gamma_idx, :) = t_values;
    dtw_d_d_val(gamma_idx, :) = d_values;
    
    % Compute GLM statistics for normalized DTW features.
    [q_values, p_values, t_values] = compute_glm_group_difference(dtw_d_norm_vec, group, covariates);
    [d_values] = compute_cohens_d_effect(dtw_d_norm_vec, group);
    dtw_d_norm_q_val(gamma_idx, :) = q_values;
    dtw_d_norm_p_val(gamma_idx, :) = p_values;
    dtw_d_norm_t_val(gamma_idx, :) = t_values;
    dtw_d_norm_d_val(gamma_idx, :) = d_values;
end

%% McNemar's sensitivity test
%% Compute Thresholded DTW Normalized Q-Values
dtw_norm_thresh = dtw_d_norm_q_val < 0.05;

% Preallocate matrices for McNemar's test p-values and sign differences.
p_val = zeros(num_gammas);
s_val = zeros(num_gammas);

% Loop over each pair of gamma values (only the lower-triangular portion is computed).
for i = 2:num_gammas
    for k = 1:(i-1)
        % Extract the binary threshold vectors for the current pair of gamma values.
        x = dtw_norm_thresh(i, :);
        y = dtw_norm_thresh(k, :);
        
        % Compute the p-value and direction (sign) using McNemar's test.
        [p, s] = compute_mcnemar_stats(x, y);
        
        % Store the results in the corresponding matrices.
        p_val(i, k) = p;
        s_val(i, k) = s;
    end
end
%% FDR Correction and Thresholding
% Vectorize the p-values, perform Benjamini-Hochberg FDR correction,
% and then convert the corrected q-values back into matrix form.
q_val = mafdr(icatb_mat2vec(p_val), "BHFDR", true);

% Compute a thresholded statistic: -log10(q) multiplied by the sign value.
q_val_thresh = icatb_vec2mat(-log10(q_val) .* icatb_mat2vec(s_val));

% Also create a binary matrix indicating whether the corrected q-value is below 0.05.
q_val_threshold = icatb_vec2mat(q_val < 0.05);

%% 

