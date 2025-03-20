%% Clear Environment and Set Paths
clc;clear;
addpath(genpath('./'));

%% Data Initialization
% Define parameters for simulated fMRI time courses.
num_subjects  = 100;         % Number of subjects
num_timepoints = 150;         % Number of timepoints per subject
num_components = 10;         % Number of components per subject
num_features   = nchoosek(num_components, 2); % Number of unique component pairs
num_sessions   = 2;          % Number of sessions
Tr             = 0.72;       % Repetition time (s)
Fs             = 1/Tr;       % Sampling frequency (Hz)
nyquist        = Fs/2;       % Nyquist frequency
band           = [0.01 0.198];  % Frequency band for post-processing
cutoff         = [0.6 1.4];     % Filter cutoff parameters

% Generate example post-processed fMRI timecourses for each session.
% The post_processing function is assumed to be available.
subTcs_sess_1 = randn(num_subjects, num_timepoints, num_components);
subTcs_sess_1 = post_processing(subTcs_sess_1, Tr, band, cutoff);

subTcs_sess_2 = randn(num_subjects, num_timepoints, num_components);
subTcs_sess_2 = post_processing(subTcs_sess_2, Tr, band, cutoff);

% Combine both sessions into a single 4D matrix with dimensions:
% [num_sessions x num_subjects x num_timepoints x num_components]
subTcs_sess_1_reshaped = reshape(subTcs_sess_1, [1, size(subTcs_sess_1)]);
subTcs_sess_2_reshaped = reshape(subTcs_sess_2, [1, size(subTcs_sess_2)]);
subTcs = cat(1, subTcs_sess_1_reshaped, subTcs_sess_2_reshaped);

%% DTW Parameters
% Define a range of gamma values for the DTW computation.
gamma_vec = 0.25:0.25:3;
num_gammas = length(gamma_vec);

%% DTW Computation Loop
% This loop computes the DTW distance and normalized DTW distance for every 
% unique component pair across all subjects, sessions, and gamma values. 
% The results are stored in 4D arrays: 
%   full_dtw_d       -> [num_gammas x num_sessions x num_subjects x num_features]
%   full_dtw_d_norm  -> [num_gammas x num_sessions x num_subjects x num_features]
%
% "icatb_mat2vec" is used to vectorize the component-wise DTW matrix into a vector
% corresponding to the unique component pairs.

% Preallocate output matrices
full_dtw_d = zeros(num_gammas, num_sessions, num_subjects, num_features);
full_dtw_d_norm = zeros(num_gammas, num_sessions, num_subjects, num_features);

for gamma_idx = 1:num_gammas
    gamma = gamma_vec(gamma_idx);
    
    for sess = 1:num_sessions
        for sub = 1:num_subjects
            % Display progress information in the command window.
            fprintf('Processing gamma = %.2f(%d/%d), session = %d/%d, subject = %d/%d\n', ...
                gamma, gamma_idx, num_gammas, sess, num_sessions, sub, num_subjects);
            
            % Initialize matrices to store DTW distances for the current subject.
            % These matrices are of size [num_components x num_components] where
            % only the lower triangular part (unique component pairs) is filled.
            sub_dtw_d = zeros(num_components, num_components);
            sub_dtw_d_norm = zeros(num_components, num_components);
            
            % Loop over all unique component pairs (comp1 > comp2)
            for comp1 = 2:num_components
                for comp2 = 1:(comp1 - 1)
                    % Extract the time courses for the two components and z-score them.
                    x = squeeze(subTcs(sess, sub, :, comp1));
                    y = squeeze(subTcs(sess, sub, :, comp2));
                    x = zscore(x);
                    y = zscore(y);
                    
                    % Compute DTW metrics using the custom DTW framework.
                    [dtw_d, dtw_d_norm] = DTWFramework(x, y, gamma, Tr, band);
                    
                    % Store the computed DTW distances in the lower-triangular matrix.
                    sub_dtw_d(comp1, comp2) = dtw_d;
                    sub_dtw_d_norm(comp1, comp2) = dtw_d_norm;
                end
            end
            
            % Vectorize the lower-triangular matrices and store them in the 4D arrays.
            full_dtw_d(gamma_idx, sess, sub, :) = icatb_mat2vec(sub_dtw_d);
            full_dtw_d_norm(gamma_idx, sess, sub, :) = icatb_mat2vec(sub_dtw_d_norm);
        end
    end
end

%% Compute Test-Retest Components for DTW
% Preallocate arrays to store between-session DTW statistics for each gamma
% and each unique component pair (feature). Both raw and normalized DTW metrics 
% are computed.
dtw_d_trt       = zeros(num_gammas, num_features); % Composite TRT score (raw DTW)
dtw_d_norm_trt  = zeros(num_gammas, num_features); % Composite TRT score (normalized DTW)

dtw_d_wsv       = zeros(num_gammas, num_features); % Within-subject variability (WSV, raw DTW)
dtw_d_norm_wsv  = zeros(num_gammas, num_features); % Within-subject variability (WSV, normalized DTW)

dtw_d_bsv       = zeros(num_gammas, num_features); % Between-subject variability (BSV, raw DTW)
dtw_d_norm_bsv  = zeros(num_gammas, num_features); % Between-subject variability (BSV, normalized DTW)

% Loop over each gamma value.
for gamma_idx = 1:num_gammas
    gamma = gamma_vec(gamma_idx);
    fprintf('Processing gamma: %.2f (%d/%d)\n', gamma, gamma_idx, num_gammas);
    
    % Loop over each feature (unique component pair).
    for feat_num = 1:num_features
        % Extract DTW distances for session 1 and session 2 (raw DTW)
        sess_1_dtw_d = squeeze(full_dtw_d(gamma_idx, 1, :, feat_num));
        sess_2_dtw_d = squeeze(full_dtw_d(gamma_idx, 2, :, feat_num));
        
        % Compute composite test-retest metrics (TRT, WSV, BSV) for raw DTW.
        [trt, wsv, bsv] = compute_composite_test_retest(sess_1_dtw_d, sess_2_dtw_d);
        dtw_d_trt(gamma_idx, feat_num) = trt;
        dtw_d_wsv(gamma_idx, feat_num) = wsv;
        dtw_d_bsv(gamma_idx, feat_num) = bsv;
        
        % Extract DTW distances for session 1 and session 2 (normalized DTW)
        sess_1_dtw_d_norm = squeeze(full_dtw_d_norm(gamma_idx, 1, :, feat_num));
        sess_2_dtw_d_norm = squeeze(full_dtw_d_norm(gamma_idx, 2, :, feat_num));
        
        % Compute composite test-retest metrics (TRT, WSV, BSV) for normalized DTW.
        [trt, wsv, bsv] = compute_composite_test_retest(sess_1_dtw_d_norm, sess_2_dtw_d_norm);
        dtw_d_norm_trt(gamma_idx, feat_num) = trt;
        dtw_d_norm_wsv(gamma_idx, feat_num) = wsv;
        dtw_d_norm_bsv(gamma_idx, feat_num) = bsv;
    end
end

%% 

function [trt, wsv, bsv] = compute_composite_test_retest(sess_1, sess_2)
% compute_composite_test_retest - Computes a composite test-retest metric.
%
% Syntax:
%   [trt, wsv, bsv] = compute_composite_test_retest(sess_1, sess_2)
%
% Inputs:
%   sess_1 - Vector of measurements for session 1.
%   sess_2 - Vector of measurements for session 2.
%
% Outputs:
%   trt    - Composite test-retest metric computed as 1 - (wsv/bsv).
%   wsv    - Parametric statistic (weighted statistic) given by the absolute t-statistic
%            from a paired t-test between sess_1 and sess_2.
%   bsv    - Nonparametric statistic (BST) computed from the median absolute deviation
%            scaled by the interquartile range.
%
% Description:
%   The function performs a paired t-test between two sessions to obtain the t-statistic,
%   whose absolute value is taken as the weighted statistic (WSV). In addition, it computes 
%   a nonparametric statistic (BST) based on the absolute deviations of the measurements 
%   from their means. The composite test-retest metric (TRT) is then defined as:
%
%       trt = 1 - (wsv / bsv)
%
% A smaller ratio (i.e., a lower WSV relative to BST) leads to a higher TRT value, 
% indicating better test-retest reliability.
%
% Example:
%   % Suppose sess1 and sess2 are vectors of measurements:
%   [trt, wsv, bsv] = compute_composite_test_retest(sess1, sess2);
%

    % Perform a paired t-test between the two sessions.
    [~, ~, ~, stats] = ttest(sess_1, sess_2);
    wsv = abs(stats.tstat);
    
    % Compute nonparametric deviations from the mean for each session.
    nonparam_sess_1 = abs(sess_1 - mean(sess_1));
    nonparam_sess_2 = abs(sess_2 - mean(sess_2));
    
    % Combine deviations from both sessions into a single vector.
    nonparam = [nonparam_sess_1(:); nonparam_sess_2(:)];
    
    % Compute the BST: median divided by the scaled interquartile range.
    bsv = median(nonparam) / (iqr(nonparam) / sqrt(length(nonparam)));
    
    % Compute the composite test-retest metric.
    trt = 1 - (wsv / bsv);
end
