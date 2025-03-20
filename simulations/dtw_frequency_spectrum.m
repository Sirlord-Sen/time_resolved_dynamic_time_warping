%% DTW and Correlation Analysis on Filtered White Noise
% This script investigates the effect of low-pass filtering on white noise
% signals by computing similarity measures between pairs of filtered signals.
%
% The analysis includes:
%   - Generating white noise signals and detrending them.
%   - Designing low-pass filters with cutoff frequencies from 0.01 to 0.66 Hz.
%   - Computing DTW distances, normalized DTW distances, and Pearson correlations
%     for multiple realizations.
%   - Averaging the similarity measures across realizations and plotting them.
%   - Demonstrating signal alignment for a chosen filter realization.
%
% Author: Sir-Lord

%% Parameters and Initialization
num_time_points = 1000;          % Number of time points in each signal
sampling_time = 0.72;            % Sampling time (sec)
filters = 0.01:0.01:0.66;         % Array of low-pass filter cutoff frequencies (Hz)
num_filters = length(filters);   % Number of filter cutoff values
num_realizations = 100;          % Number of realizations for the main analysis

% Preallocate matrices to store DTW and correlation values
dtw_matrix = zeros(num_realizations, num_filters);
dtw_norm_matrix = zeros(num_realizations, num_filters);
corr_matrix = zeros(num_realizations, num_filters);

% Sampling frequency
Fs = 1/sampling_time; 

%% Compute Upper-Bound DTW Distances at Nyquist
% These upper-bound values are computed using an extreme filter cutoff (0.67 Hz)
% over a large number of realizations. They are later used to normalize DTW distances.
num_realizations_upper = 1000; 
dtw_norm_upper = zeros(num_realizations_upper, 1);
dtw_upper = zeros(num_realizations_upper, 1);

for i = 1:num_realizations_upper
    % Generate two independent white noise signals and remove linear trend
    signal1 = detrend(randn(num_time_points, 1));
    signal2 = detrend(randn(num_time_points, 1));
    
    % Set an extreme cutoff frequency (0.67 Hz) for upper-bound estimation
    cutoff = 0.67;
    % Determine filter order and normalized cutoff using buttord
    [n, Wn] = buttord(cutoff/(Fs/2), (cutoff+0.02)/(Fs/2), 3, 30);
    [b, a] = butter(n, Wn, 'low');  % Design the low-pass filter
    
    % Apply the filter to both signals
    filtered_signal1 = filter(b, a, signal1);
    filtered_signal2 = filter(b, a, signal2);
    
    % Z-score the filtered signals
    filtered_signal1_zscored = zscore(filtered_signal1);
    filtered_signal2_zscored = zscore(filtered_signal2);
    
    % Compute DTW distance using MATLAB's built-in dtw function
    [dist, ix] = dtw(filtered_signal1_zscored, filtered_signal2_zscored);
    dtw_upper(i) = dist;
    dtw_norm_upper(i) = dist / length(ix);
end

%% Main Loop: Compute Similarity Metrics Across Realizations and Filters
% Preallocate storage for DTW distances and store filtered signal pairs.
dtw_distances = zeros(num_realizations, num_filters);
signal_pairs = zeros(num_filters, 2, num_time_points);

for i = 1:num_realizations
    fprintf('Realization %d/%d\n', i, num_realizations);
    
    % Generate two new white noise signals and detrend them
    signal1 = detrend(randn(num_time_points, 1));
    signal2 = detrend(randn(num_time_points, 1));
    
    % Loop over each low-pass filter cutoff value
    for j = 1:num_filters
        cutoff = filters(j);
        
        % Design low-pass filter using buttord for current cutoff
        [n, Wn] = buttord(cutoff/(Fs/2), (cutoff+0.02)/(Fs/2), 3, 30);
        [b, a] = butter(n, Wn, 'low');
        
        % Apply filter to both signals
        filtered_signal1 = filter(b, a, signal1);
        filtered_signal2 = filter(b, a, signal2);
        
        % Z-score the filtered signals
        filtered_signal1_zscored = zscore(filtered_signal1);
        filtered_signal2_zscored = zscore(filtered_signal2);
        
        % Store the filtered signal pair for later visualization
        signal_pairs(j, :, :) = [filtered_signal1_zscored, filtered_signal2_zscored]';
        
        % Compute DTW distance between the two filtered signals
        [dist, ix] = dtw(filtered_signal1_zscored, filtered_signal2_zscored);
        dtw_distances(i, j) = dist;
        
        % Normalize DTW distances using the upper-bound values
        dtw_norm_matrix(i, j) = 1 - (dist/length(ix)) / mean(dtw_norm_upper);
        dtw_matrix(i, j) = 1 - dist / mean(dtw_upper);
        
        % Compute Pearson correlation between the two filtered signals
        corr_matrix(i, j) = corr(filtered_signal1_zscored, filtered_signal2_zscored);
    end
end

% Average the similarity measures over all realizations
avg_dtw = mean(dtw_matrix, 1);        % Average normalized DTW (based on raw distance)
avg_dtw_norm = mean(dtw_norm_matrix, 1);% Average normalized DTW (based on path length)
avg_corr = mean(corr_matrix, 1);        % Average correlation

%% Plot Average Similarity Measures vs. Filter Cutoff
figure;
hold on;
plot(filters, avg_corr, '-o', 'Color', [0 148 17]/255, 'LineWidth', 7, 'DisplayName', 'Correlation');
plot(filters, avg_dtw, '-o', 'Color', [1 0 0], 'LineWidth', 7, 'DisplayName', 'DTW');
plot(filters, avg_dtw_norm, '-o', 'Color', [0 0 1], 'LineWidth', 7, 'DisplayName', 'Normalized DTW');
hold off;
xlim([0 0.67]);
set(gca, 'FontSize', 35, 'LineWidth', 1.5);
grid on;
legend('Orientation','horizontal','Location','best');
fig = gcf;
fig.Position = [100, 100, 1600, 900];


%% Plot Aligned Signals for a Specific Filter Realization
% Choose a filter index (e.g., f_idx = 10) to visualize signal alignment.
f_idx = 10;
x_aligned = squeeze(signal_pairs(f_idx, 1, :));
y_aligned = squeeze(signal_pairs(f_idx, 2, :));
[dtw_d, ix, iy] = dtw(x_aligned, y_aligned);

% Plot the aligned signals using the DTW mapping indices.
hF = figure();
hA = axes(hF);
plot(x_aligned(ix), 'Color', [98 75 131]/255, 'LineWidth', 7);
hold on;
plot(y_aligned(iy), 'Color', [158 151 36]/255, 'LineWidth', 7);
hold off;
xlim([0 length(ix)]);
box off;
grid off;
xticks([]);
yticks([]);
set(hA, 'XColor', 'none', 'YColor', 'none');
fig = gcf;
fig.Position = [100, 100, 1600, 700];

%% Plot Raw Filtered Signals for a Specific Filter Realization
% Display the two filtered signals (without DTW alignment) for comparison.
x_raw = squeeze(signal_pairs(f_idx, 1, :));
y_raw = squeeze(signal_pairs(f_idx, 2, :));
hF = figure();
hA = axes(hF);
plot(x_raw, 'Color', [98 75 131]/255, 'LineWidth', 7);
hold on;
plot(y_raw, 'Color', [158 151 36]/255, 'LineWidth', 7);
hold off;
xlim([0 length(x_raw)]);
box off;
grid off;
xticks([]);
yticks([]);
set(hA, 'XColor', 'none', 'YColor', 'none');
fig = gcf;
fig.Position = [100, 100, 1600, 700];
