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

%% Generate Surrogate fMRI Time Courses
num_surrogates = 10;

% Generate surrogate time courses for the f1 pipeline using a null model.
surr_subTcs_f1 = pr_null_model(subTcs_f1, num_surrogates);

% Preallocate the array for the full surrogate time courses.
surr_subTcs_full = zeros(size(surr_subTcs_f1));

% Loop over each surrogate.
for surr = 1:num_surrogates
    % Process the surrogate f1 time courses with the f2 parameters to extract
    % the low-frequency component.
    surr_tmp_f1_low = post_processing(squeeze(surr_subTcs_f1(:, :, :, surr)), Fs, band_f2, cutoff_f2);
    
    % Compute the high-frequency component by subtracting the low-frequency
    % component from the original f1 data.
    surr_tmp_f1_high = subTcs_f1 - surr_tmp_f1_low;
    
    % Combine the high-frequency component from f1 with the original f2 data
    % to form the full surrogate dataset.
    surr_subTcs_full(:, :, :, surr) = subTcs_f2 + surr_tmp_f1_high;
end

%% 
% Perform DTW on the surrogate data to generate null
% Perform a wilcoxin sign rank test on the DTW metric vs null

%% Steps:
% 1. Perform DTW on the Surrogate Data to Generate Null Distribution
% 2. Perform DTW on the original signal F1
% 3. Perform a wilcoxin signrank test on output of step 2 vs step 1