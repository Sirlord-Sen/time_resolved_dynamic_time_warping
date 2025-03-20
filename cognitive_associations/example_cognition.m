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
dtw_d_norm = squeeze(full_dtw_d_norm(1, :, :));

%% Setting up Covariates
% Define group labels: schizophrenia subjects (1) and controls (0).
sz_label = 1;
cn_label = 0;
group = randi([cn_label, sz_label], num_subjects, 1);  % Random binary group assignment

% Create covariate vectors.
age  = randi([20, 80], num_subjects, 1);  % Random ages between 20 and 80
sex  = randi([0, 1], num_subjects, 1);     % Random binary sex assignment
site = randi([1, 3], num_subjects, 1);     % Random site assignment (3 possible sites)
mfd  = randn(num_subjects, 1);             % Random mean frame displacement

% Dummy encoding for site (exclude the first column to avoid collinearity).
site_dummy = dummyvar(categorical(site));
% For full covariates, include group
covariates = [group, mfd, age, sex, site_dummy(:, 2:end)];

%% Set up Covariates for Schizophrenia only
dtw_d_norm_sz = dtw_d_norm(group == sz_label, :);
age_sz   = age(group == sz_label);
sex_sz   = sex(group == sz_label);
site_sz  = site(group == sz_label);
mfd_sz   = mfd(group == sz_label);

% Dummy encoding for site (exclude the first column to avoid collinearity).
site_dummy_sz = dummyvar(categorical(site_sz));
covariates_sz = [mfd_sz, age_sz, sex_sz, site_dummy_sz(:, 2:end)];
 
%% Set up Covariates for Controls only
dtw_d_norm_cn = dtw_d_norm(group == cn_label, :);
age_cn   = age(group == cn_label);
sex_cn   = sex(group == cn_label);
site_cn  = site(group == cn_label);
mfd_cn   = mfd(group == cn_label);

% Dummy encoding for site (exclude the first column to avoid collinearity).
site_dummy_cn = dummyvar(categorical(site_cn));
covariates_cn = [mfd_cn, age_cn, sex_cn, site_dummy_cn(:, 2:end)];

%% Create PANSS Score Distribution
% Generate a normal distribution to mimic cognitive scores
cog_score = randn(num_subjects, 1);

cog_score_sz = cog_score(group == sz_label);
cog_score_cn = cog_score(group == cn_label);
%% GLM Association Analysis
% Specify the distribution as "normal" or remove to choose default "normal".
dist = "normal";

%% Perform GLM on all subjects
[q_values, p_values, t_values, b_values] = compute_glm_associations(dtw_d_norm, covariates, cog_score, dist);

%% Perform GLM on only schizophrenia
[q_values_sz, p_values_sz, t_values_sz, b_values_sz] = compute_glm_associations(dtw_d_norm_sz, covariates_sz, cog_score_sz, dist);
 
%% Perform GLM on only controls
[q_values_cn, p_values_cn, t_values_cn, b_values_cn] = compute_glm_associations(dtw_d_norm_cn, covariates_cn, cog_score_cn, dist);

