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
gamma_vec = 1;
num_gammas = length(gamma_vec);

%% DTW and nDTW Computation across all subjects and gammas.
[full_dtw_d, full_dtw_d_norm] = compute_subject_dtw_framework(subTcs, gamma_vec, Tr, band);

% For illustration, pick one gamma value (first one).
dtw_d = squeeze(full_dtw_d(1, :, :));

%% Setting Up Group Difference Covariates
% Define group labels: schizophrenia subjects (1) and controls (0).
schizophrenia_label = 1;
control_label       = 0;
group = randi([control_label, schizophrenia_label], num_subjects, 1);  % Random binary group assignment

% Create covariate vectors.
age  = randi([20, 80], num_subjects, 1);  % Random ages between 20 and 80
sex  = randi([0, 1], num_subjects, 1);     % Random binary sex assignment
site = randi([1, 3], num_subjects, 1);     % Random site assignment (3 possible sites)
mfd  = randn(num_subjects, 1);             % Random mean frame displacement

% For symptom analysis, use only schizophrenia subjects.
dtw_d_sz = dtw_d(group == schizophrenia_label, :);
age_sz   = age(group == schizophrenia_label);
sex_sz   = sex(group == schizophrenia_label);
site_sz  = site(group == schizophrenia_label);
mfd_sz   = mfd(group == schizophrenia_label);

% Dummy encoding for site (exclude the first column to avoid collinearity).
site_dummy_sz = dummyvar(categorical(site_sz));
covariates = [mfd_sz, age_sz, sex_sz, site_dummy_sz(:, 2:end)];
%% Create PANSS Score Distribution
% Generate a non-negative whole number distribution to mimic PANSS scores for schizophrenia subjects.
num_sz = sum(group == schizophrenia_label);
panss_score = randi([1, 45], num_sz, 1);
%% GLM Association Analysis
% Specify the distribution. Since PANSS scores are Poisson distributed, use "poisson".
dist = "poisson";
[q_values, p_values, t_values, b_values] = compute_glm_associations(dtw_d_sz, covariates, panss_score, dist);
