%% Script: phase_amplitude_non_stationarity.m
% This script demonstrates how to use Dynamic Time Warping (DTW) to analyze
% the effects of phase and amplitude perturbations on a sine wave signal.
%
% The script performs the following tasks:
%   1. Defines a sine wave signal and simulation parameters.
%   2. Computes DTW upper bounds for random phase and amplitude perturbations.
%   3. Simulates time-varying phase perturbations and computes DTW/correlation measures.
%   4. Simulates time-varying amplitude perturbations and computes DTW/correlation measures.
%   5. Visualizes the results and example signals.
%
% Required functions:
%   - dtw_custom.m
%   - get_dtw_phase_upper_bound.m
%   - get_dtw_amplitude_upper_bound.m
%
% Author: Sir-Lord Wiafe

%% Initialize Environment
clear; close all; clc;

%% Add Paths (update these as needed)
addpath(genpath('./'));
addpath(genpath('./utils/'));

%% Signal and Simulation Parameters
Ts = 0.72;                    % Sampling interval (sec)
Fs = 1/Ts;                    % Sampling frequency (Hz)
nyquist = Fs/2;               % Nyquist frequency
L = 1000;                     % Number of samples in the full time series
t = (0:1/Fs:(Ts*L - 1/Fs))';   % Time vector

% Sine wave parameters
f = 0.01;         % Frequency of the sine wave (Hz)
A = 1;            % Amplitude of the sine wave
phi = 0;          % Initial phase (radians)

% Perturbation parameters
amp_max = 50;     % Maximum amplitude deviation
phi_max = 2*pi;   % Maximum phase deviation (radians)

% DTW parameters
win_size = 300;   % Warping window size constraint
gamma = 1;        % Exponent for the DTW distance measure

% Generate the original sine wave signal
x = sin(2*pi*f*t + phi);

%% Compute DTW Upper Bounds for Phase Perturbations
num_trials = 50;  % Number of trials for random phase perturbations
[dtw_phase_max_dist, dtw_norm_phase_max_dist] = get_dtw_phase_upper_bound(x, win_size, gamma, f, Fs, phi_max, num_trials);
dtw_phase_upper = mean(dtw_phase_max_dist);
dtw_norm_phase_upper = mean(dtw_norm_phase_max_dist);

%% Compute DTW Upper Bounds for Amplitude Perturbations
num_trials = 50;  % Number of trials for random amplitude perturbations
[dtw_amp_max_dist, dtw_amp_norm_max_dist] = get_dtw_amplitude_upper_bound(x, win_size, gamma, f, Fs, amp_max, num_trials);
dtw_amp_upper = mean(dtw_amp_max_dist);
dtw_norm_amp_upper = mean(dtw_amp_norm_max_dist);

%% Simulation: Time-Varying Phase Perturbation
scale_phase = linspace(phi, phi_max, num_steps);
scale_phase_labels = arrayfun(@(s) num2str(floor(s * 100) / 100), scale_phase, 'UniformOutput', false);

% Pre-allocate arrays for correlation and DTW metrics
corr_sim_phase = zeros(1, num_steps);
dtw_sim_phase = zeros(1, num_steps);
dtw_norm_sim_phase = zeros(1, num_steps);

% Loop over different phase scales to simulate time-varying phase
for k = 1:num_steps
    % Create a time-varying phase: linear progression from phi to scale_phase(k)
    phi_t = phi + (scale_phase(k) - phi) * t / t(end);
    y = sin(2*pi*f*t + phi_t);
    
    % Compute correlation and DTW measures between the original and perturbed signals
    corr_sim_phase(k) = corr(zscore(x), zscore(y));
    [~, dtw_sim_phase(k), ix] = dtw_custom(zscore(x), zscore(y), win_size, gamma);
    dtw_norm_sim_phase(k) = dtw_sim_phase(k) / length(ix);
end

% Normalize DTW measures relative to the DTW phase upper bounds
dtw_sim_phase = 1 - dtw_sim_phase ./ dtw_phase_upper;
dtw_norm_sim_phase = 1 - dtw_norm_sim_phase ./ dtw_norm_phase_upper;

%% Plot Results: Time-Varying Phase Perturbation
figure(1);
plot(corr_sim_phase, 'LineWidth', 7, 'Color', [0 148 17]/255);
hold on;
plot(dtw_sim_phase, 'LineWidth', 7, 'Color', [1 0 0]);
plot(dtw_norm_sim_phase, 'LineWidth', 7, 'Color', [0 0 1]);
hold off;
set(gca, 'FontSize', 35, 'LineWidth', 3);
xticks(1:(num_steps+1));
xticklabels({scale_phase_labels{:}});
legend('Correlation', 'DTW', 'nDTW', 'Location','best')
grid on;
box off;
% Adjust figure size (optional)
fig = gcf;
fig.Position = [100, 100, 1600, 700];

%% Visualize Example Signals for Phase Perturbation
hF = figure();
hA = axes(hF);
plot(zscore(x), 'LineWidth', 10, 'Color', [98 75 131]/255);
hold on;
plot(zscore(y), 'LineWidth', 10, 'Color', [158 151 36]/255);
hold off;
box off;
grid off;
xticks([]);
yticks([]);
set(hA, 'XColor', 'none', 'YColor', 'none');
fig = gcf;
fig.Position = [100, 100, 1600, 700];

%% Simulation: Time-Varying Amplitude Perturbation
num_steps = 10;
% Combining two amplitude scales for broader range visualization
scale_amp = [linspace(A, amp_max, num_steps)];
scale_amp_labels = arrayfun(@(s) num2str(floor(s * 10) / 100), scale_amp, 'UniformOutput', false);

% Pre-allocate arrays for correlation and DTW metrics for amplitude perturbations
corr_sim_amp = zeros(1, length(scale_amp));
dtw_sim_amp = zeros(1, length(scale_amp));
dtw_norm_sim_amp = zeros(1, length(scale_amp));

% Loop over different amplitude scales to simulate time-varying amplitude
for k = 1:length(scale_amp)
    % Create a time-varying amplitude: linear modulation from A to scale_amp(k)
    A_t = A + (scale_amp(k) - A) * t / t(end);
    y = A_t .* x;
    
    % Compute correlation and DTW measures between the original and amplitude-modulated signals
    corr_sim_amp(k) = corr(zscore(x), zscore(y));
    [~, dtw_sim_amp(k), ix] = dtw_custom(zscore(x), zscore(y), win_size, gamma);
    dtw_norm_sim_amp(k) = dtw_sim_amp(k) / (length(ix) - 1);
end

% Normalize DTW measures relative to the DTW amplitude upper bounds
dtw_sim_amp = 1 - dtw_sim_amp ./ dtw_amp_upper;
dtw_norm_sim_amp = 1 - dtw_norm_sim_amp ./ dtw_norm_amp_upper;

%% Plot Results: Time-Varying Amplitude Perturbation
figure;
plot(corr_sim_amp, 'LineWidth', 7, 'Color', [0 148 17]/255);
hold on;
plot(dtw_sim_amp, 'LineWidth', 7, 'Color', [1 0 0]);
plot(dtw_norm_sim_amp, 'LineWidth', 7, 'Color', [0 0 1]);
hold off;
set(gca, 'FontSize', 35, 'LineWidth', 3);
xticks(1:(length(scale_amp)+1));
xticklabels({scale_amp_labels{:}});
ylim([0 1]);
legend('Correlation', 'DTW', 'nDTW', 'Location','best')
grid on;
box off;
fig = gcf;
fig.Position = [100, 100, 1600, 700];

%% Visualize Example Signals for Amplitude Perturbation
hF = figure();
hA = axes(hF);
plot(zscore(x), 'LineWidth', 10, 'Color', [98 75 131]/255);
hold on;
plot(zscore(y), 'LineWidth', 10, 'Color', [158 151 36]/255);
hold off;
box off;
grid off;
xticks([]);
yticks([]);
set(hA, 'XColor', 'none', 'YColor', 'none');
fig = gcf;
fig.Position = [100, 100, 1600, 700];
