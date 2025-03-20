%% Phase Modulation Analysis Using DTW
% This script evaluates how phase perturbations applied to a base sine wave
% affect various similarity indices. The analysis is performed over a range
% of modulation scales and for different DTW gamma values.
%
% The following similarity indices are computed:
%   - Correlation (via corr2)
%   - Coherence (via mscohere)
%   - DTW similarity index (1 - normalized DTW distance)
%   - Normalized DTW similarity index
%
% DTW upper bounds are computed using an extremely perturbed phase version
% of the signal. The results are plotted as functions of the modulation scale,
% with the DTW values for gamma = 1 highlighted in the first figure.
%
% Author: Sir-Lord

%% Setup Base Signal and Parameters
Tr = 0.72;                   % Sampling interval (s)
t = 0:Tr:1500;               % Time vector
f = 0.01;                    % Base frequency (Hz)
fm = 0.0001;                 % Modulation frequency (Hz)

% Create a base sine wave
signal1 = sin(2*pi*f*t);

% Define DTW gamma values and modulation scales (log-spaced)
gamma_vec = [0.25; 0.5; 1; 1.5; 2];
scale_vec = logspace(log10(0.5), log10(15), 20);

% Preallocate matrices for storing similarity indices
dtw_scale     = zeros(length(scale_vec), length(gamma_vec));  % DTW similarity index
dtw_norm_scale = zeros(length(scale_vec), length(gamma_vec)); % Normalized DTW similarity index
corr_scale    = zeros(length(scale_vec), 1);                  % Correlation similarity index
coh_scale     = zeros(length(scale_vec), 1);                  % Coherence similarity index

%% Compute Similarity Indices for Phase Modulation
% For each gamma value, first compute an upper-bound DTW distance using an
% extremely perturbed phase signal. Then, loop over modulation scales to
% calculate the similarity indices.
for gamma_idx = 1:length(gamma_vec)
    fprintf('Processing gamma index %d/%d\n', gamma_idx, length(gamma_vec));
    gamma = gamma_vec(gamma_idx);
    
    % Compute DTW upper bound using an extremely perturbed phase signal
    signal1_z = zscore(signal1);
    % Use a very large phase perturbation for upper bound estimation
    signal2_z = zscore(sin(2*pi*f*t + 100000000*sin(2*pi*fm*t)));
    [~, upper_dist, ix] = dtw_custom(signal1_z, signal2_z, 1000, gamma);
    upper_norm_dist = upper_dist / length(ix);
    
    % Loop over modulation scales
    for k = 1:length(scale_vec)
        % Apply phase modulation: modulation is a slow sine wave scaled by scale_vec(k)
        modulation = scale_vec(k) * sin(2*pi*fm*t);
        signal2 = sin(2*pi*f*t + modulation);
        
        % Z-score the base and modulated signals
        signal1_z = zscore(signal1);
        signal2_z = zscore(signal2);
        
        % Compute correlation similarity (using 2-D correlation)
        corr_scale(k) = corr2(signal1_z, signal2_z);
        
        % Compute coherence similarity (average coherence over frequencies)
        tmp_coh = mscohere(signal1_z, signal2_z);
        coh_scale(k) = mean(tmp_coh);
        
        % Compute DTW distance between z-scored signals
        [~, dist, ix] = dtw_custom(signal1_z, signal2_z, 1000, gamma);
        
        % Normalize DTW distances relative to the upper bound
        dtw_scale(k, gamma_idx) = 1 - dist / upper_dist;
        dtw_norm_scale(k, gamma_idx) = 1 - (dist / length(ix)) / upper_norm_dist;
    end
end

%% Plot Similarity Indices (Gamma = 1 Highlighted)
figure;
hold on;
plot(scale_vec, corr_scale, '-o', 'Color', [0 148 17]/255, 'LineWidth', 7, 'DisplayName', 'Correlation');
plot(scale_vec, coh_scale, '-o', 'Color', [246 154 25]/255, 'LineWidth', 7, 'DisplayName', 'Coherence');
% Plot DTW and normalized DTW for gamma = 1 (third entry in gamma_vec)
plot(scale_vec, dtw_scale(:, 3), '-o', 'Color', [1 0 0], 'LineWidth', 7, 'DisplayName', 'DTW (\gamma=1)');
plot(scale_vec, dtw_norm_scale(:, 3), '-o', 'Color', [0 0 1], 'LineWidth', 7, 'DisplayName', 'nDTW (\gamma=1)');
hold off;

set(gca, 'XScale', 'log');   % Use logarithmic x-axis for modulation scale
set(gca, 'FontSize', 35, 'LineWidth', 1.5);
grid on;
legend('Orientation','horizontal','Location','best');
fig = gcf;
fig.Position = [100, 100, 1600, 700];


%% Compute and Plot Rate of Change in DTW Similarity
% Define colors for each gamma value
colors_vec = [ ...
    0.3010 0.7450 0.9330; ...
    0 0.4470 0.7410; ...
    1.0 0.0 0.0; ...
    0.8 0.6 0.4;  ...
    0.5 0.3 0.2];

% Calculate absolute rate of change (first difference) along the modulation scale
dtw_rate_of_change = abs(diff(dtw_scale, 1, 1));
% Pad with a row of zeros to match the original size
dtw_rate_of_change = [dtw_rate_of_change; zeros(1, size(dtw_rate_of_change, 2))];

dtw_norm_rate_of_change = abs(diff(dtw_norm_scale, 1, 1));
dtw_norm_rate_of_change = [dtw_norm_rate_of_change; zeros(1, size(dtw_norm_rate_of_change, 2))];

figure;
hold on;
% Plot rate of change for each gamma value
for gamma_idx = 1:length(gamma_vec)
    % Plot DTW rate of change (solid line)
    plot(scale_vec, dtw_rate_of_change(:, gamma_idx), '-o', 'LineWidth', 7, ...
         'Color', colors_vec(gamma_idx, :), ...
         'DisplayName', ['DTW \gamma=' num2str(gamma_vec(gamma_idx))]);
    
    % Plot normalized DTW rate of change (dashed line)
    plot(scale_vec, dtw_norm_rate_of_change(:, gamma_idx), ':o', 'LineWidth', 7, ...
         'Color', colors_vec(gamma_idx, :), ...
         'DisplayName', ['nDTW \gamma=' num2str(gamma_vec(gamma_idx))]);
end

set(gca, 'XScale', 'log');   % Log scale for modulation scale
set(gca, 'YScale', 'log');   % Log scale for y-axis to emphasize small differences
set(gca, 'FontSize', 35, 'LineWidth', 1.5);
grid on;
hold off;

% Manually create custom legend entries in a 2x5 layout (2 rows, 5 columns)
legend_entries = cell(2 * length(gamma_vec), 1);
for gamma_idx = 1:length(gamma_vec)
    legend_entries{2*gamma_idx - 1} = ['DTW \gamma=' num2str(gamma_vec(gamma_idx))];
    legend_entries{2*gamma_idx}     = ['nDTW \gamma=' num2str(gamma_vec(gamma_idx))];
end

lgd = legend(legend_entries);
lgd.FontSize = 20;
lgd.NumColumns = 5;
lgd.Location = 'north';

fig = gcf;
fig.Position = [100, 100, 1600, 900];
