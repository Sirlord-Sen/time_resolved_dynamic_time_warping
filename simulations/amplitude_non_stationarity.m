%% DTW Sensitivity Analysis for Amplitude Non-Stationarity
% This script evaluates how different similarity indices respond to slow 
% amplitude modulations applied to a base sine wave. Specifically, it computes:
%   - Correlation (using 2-D correlation, via corr2)
%   - Coherence (via mscohere)
%   - DTW and Normalized DTW similarity indices
%
% DTW distances are computed using varying gamma exponents (which weight the 
% local distance measure) and are normalized relative to an upper-bound computed 
% from an extremely perturbed version of the signal.
%
% Required external functions:
%   - dtw_custom.m
%
%
% Author: Sir-Lord

%% Initialize Environment and Define Paths
clear; close all; clc;
addpath(genpath('./'));
addpath(genpath('./utils/'));
%% Create Base Sine Wave Signal
Tr = 0.72;              % Sampling interval (s)
t = 0:Tr:1500;          % Time vector (s)
f = 0.01;               % Base frequency (Hz)
fm = 0.0001;            % Modulation frequency (Hz)
signal1 = sin(2*pi*f*t);% Base sine wave signal (frequency = 0.01 Hz)

%% Define Analysis Parameters
gamma_vec = [0.25; 0.5; 1; 1.5; 2];  % DTW gamma values (exponent for local distance)
scale_vec = logspace(log10(0.01), log10(100), 20); % Modulation scales (logarithmic)

% Pre-allocate arrays for storing similarity measures as functions of modulation scale.
% For DTW measures, we store results for each gamma value.
dtw_scale     = zeros(length(scale_vec), length(gamma_vec));  % DTW similarity index
dtw_norm_scale = zeros(length(scale_vec), length(gamma_vec));  % Normalized DTW similarity index

% Correlation and coherence are computed once per modulation scale (the last gamma loop's results)
corr_scale = zeros(length(scale_vec), 1);  % Correlation similarity
coh_scale  = zeros(length(scale_vec), 1);    % Coherence similarity

%% Compute Upper-Bound DTW Distance for Extreme Perturbation
% For each gamma value, compute an upper bound on the DTW distance using an extremely 
% perturbed version of the signal. This bound will be used to normalize DTW distances.
for gamma_idx = 1:length(gamma_vec)
    fprintf('Processing gamma index: %d\n', gamma_idx);
    gamma = gamma_vec(gamma_idx);
    
    % Compute z-scored base signal
    signal1_z = zscore(signal1);
    
    % Generate an extremely perturbed signal using a very large amplitude perturbation.
    extreme_scale = 1 + 1e10 * sin(2*pi*fm*t);
    signal2_z = zscore(signal1 .* extreme_scale);
    
    % Compute DTW distance between the base and extremely perturbed signals
    [~, upper_dist, ix] = dtw_custom(signal1_z, signal2_z, [], gamma);
    upper_norm_dist = upper_dist / length(ix);
    
    %% For Each Modulation Scale, Compute Similarity Indices
    % Loop over a set of amplitude modulation scales to create a time-varying amplitude.
    for k = 1:length(scale_vec)
        % Create slow amplitude modulation: modulation varies linearly with sin(2*pi*fm*t)
        modulation = 1 + scale_vec(k) * sin(2*pi*fm*t);
        signal2 = signal1 .* modulation;  % Apply modulation to the base signal
        
        % Z-score both signals for robust comparison
        signal1_z = zscore(signal1);
        signal2_z = zscore(signal2);
        
        % Compute correlation similarity (using 2-D correlation since signals are vectors)
        corr_scale(k) = corr2(signal1_z, signal2_z);
        
        % Compute coherence using mscohere; take the average of the coherence spectrum
        tmp_coh = mscohere(signal1_z, signal2_z);
        coh_scale(k) = mean(tmp_coh);
        
        % Compute DTW distance between the z-scored signals
        [~, dist, ix] = dtw_custom(signal1_z, signal2_z, [], gamma);
        
        % Normalize DTW distances relative to the upper bound
        dtw_scale(k, gamma_idx)     = 1 - dist / upper_dist;
        dtw_norm_scale(k, gamma_idx) = 1 - (dist / length(ix)) / upper_norm_dist;
    end
end

%% Plot Similarity Indices as a Function of Modulation Scale
% In this plot, we use the DTW and normalized DTW results for gamma = 1 (third entry in gamma_vec)
% for comparison with correlation and coherence measures.
figure;
hold on;
plot(scale_vec, corr_scale, '-o', 'Color', [0 148 17]/255, 'LineWidth', 7);
plot(scale_vec, coh_scale, '-o', 'Color', [246 154 25]/255, 'LineWidth', 7);
plot(scale_vec, dtw_scale(:, 3), '-o', 'Color', [1 0 0], 'LineWidth', 7);
plot(scale_vec, dtw_norm_scale(:, 3), '-o', 'Color', [0 0 1], 'LineWidth', 7);
hold off;

set(gca, 'XScale', 'log');  % Use logarithmic x-axis for modulation scale
set(gca, 'FontSize', 35, 'LineWidth', 1.5);
grid on;
box off;
legend('Correlation', 'Coherence', 'DTW', 'Normalized DTW', 'Orientation', 'horizontal', 'Location', 'best');

% Adjust figure size
fig = gcf;
fig.Position = [100, 100, 1600, 700];

%% 
%% Plot Rate of Change in DTW Similarity vs. Modulation Scale
% This section calculates the absolute rate of change (first difference)
% for both DTW and normalized DTW similarity indices, and then plots these
% rates on a log-log scale for different gamma values.

% Define colors for each gamma value
colors_vec = [
    0.3010 0.7450 0.9330;
    0 0.4470 0.7410;
    1.0 0.0 0.0;
    0.8 0.6 0.4;  
    0.5 0.3 0.2;
];

% Calculate the absolute rate of change for DTW similarity
dtw_rate_of_change = abs(diff(dtw_scale, 1, 1));
% Pad with a row of zeros to match the original size
dtw_rate_of_change = [dtw_rate_of_change; zeros(1, size(dtw_rate_of_change, 2))];

% Calculate the absolute rate of change for normalized DTW similarity
dtw_norm_rate_of_change = abs(diff(dtw_norm_scale, 1, 1));
% Pad with a row of zeros to match the original size
dtw_norm_rate_of_change = [dtw_norm_rate_of_change; zeros(1, size(dtw_norm_rate_of_change, 2))];

figure;
hold on;
% Loop through each gamma value to plot DTW and normalized DTW rates
for gamma_idx = 1:length(gamma_vec)
    % Plot DTW rate of change with a solid line
    plot(scale_vec, dtw_rate_of_change(:, gamma_idx), '-o', 'LineWidth', 7, ...
         'Color', colors_vec(gamma_idx, :), ...
         'DisplayName', ['DTW \gamma=' num2str(gamma_vec(gamma_idx))]);
    
    % Plot normalized DTW rate of change with a dashed line
    plot(scale_vec, dtw_norm_rate_of_change(:, gamma_idx), ':o', 'LineWidth', 7, ...
         'Color', colors_vec(gamma_idx, :), ...
         'DisplayName', ['nDTW \gamma=' num2str(gamma_vec(gamma_idx))]);
end

% Set axes to logarithmic scale to emphasize small differences
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
set(gca, 'FontSize', 35, 'LineWidth', 1.5);
ylim([0.00003 0.8]);
grid on;
hold off;

% Create custom legend entries for clarity
legend_entries = cell(2 * length(gamma_vec), 1);
for gamma_idx = 1:length(gamma_vec)
    legend_entries{2 * gamma_idx - 1} = ['DTW \gamma=' num2str(gamma_vec(gamma_idx))];
    legend_entries{2 * gamma_idx}     = ['nDTW \gamma=' num2str(gamma_vec(gamma_idx))];
end

% Create and format the legend in a 2x5 layout (2 rows, 5 columns)
lgd = legend(legend_entries);
lgd.FontSize = 20;
lgd.NumColumns = 5;
lgd.Location = 'north';

% Adjust the figure size
fig = gcf;
fig.Position = [100, 100, 1600, 900];
