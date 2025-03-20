%% DTW Alignment and Interpolation of Subcortical Time Courses
% This script loads subcortical time course data and associated labels,
% selects two specific components, and computes their dynamic time warping (DTW)
% alignment. It then interpolates the DTW mapping to generate reference vectors 
% and visualizes both the raw and warped signals.
%
% Author: Sir-Lord Wiafe

%% Initialization and Data Loading
clear; close all; clc;

% Extract and z-score the time courses for the two selected components
x = zscore(randn(157, 1));
y = zscore(randn(157, 1));
%% Set Parameters for DTW Alignment
win_size = 45;          % DTW warping window size

%% Plot Original Raw Signals
L = length(x);
hF = figure();
hA = axes(hF);
plot(x(1:80), 'LineWidth', 15, 'Color', [98 75 131]/255);
hold on;
plot(y(1:80), 'LineWidth', 15, 'Color', [157 180 124]/255);
yline(0, '--k', 'LineWidth', 5);
hold off;
ylabel('Amplitude', 'FontSize', 13);

box off; grid off;
xticks([]); yticks([]);
set(get(hA, 'XAxis'), 'Visible', 'off');
set(get(hA, 'YAxis'), 'Visible', 'off');

fig = gcf;
fig.Position = [100, 100, 1800, 700];  % Adjust figure size


%% Compute DTW Alignment and Warp the Signals
% Compute the DTW distance and obtain the optimal alignment indices
[dtw_d, ix, iy] = dtw(x, y, win_size);

% Generate a sampling vector spanning the original time course length
x_samp = linspace(1, length(ix), L);

% Interpolate the warped signals using PCHIP interpolation based on DTW mapping
x_w = interp1(x(ix), x_samp, "pchip");
y_w = interp1(y(iy), x_samp, "pchip");

% Plot the warped signals (aligned by DTW indices)
hF = figure();
hA = axes(hF);
plot(x_w(1:80), 'LineWidth', 15, 'Color', [246 154 25]/255);
hold on;
plot(y_w(1:80), 'LineWidth', 15, 'Color', [25 172 246]/255);
yline(0, '--k', 'LineWidth', 5);
hold off;
ylabel('Amplitude', 'FontSize', 13);
% Optionally, add legend or title
box off; grid off;
xticks([]); yticks([]);
set(get(hA, 'XAxis'), 'Visible', 'off');
set(get(hA, 'YAxis'), 'Visible', 'off');

fig = gcf;
fig.Position = [100, 100, 1800, 700];

%% Compute and Plot Interpolated Reference Vectors
% Use a custom function to compute reference vectors based on the DTW mapping.
% The function "interpolate_vectors" returns the interpolated reference vectors.
[x_ref_v_inc, y_ref_v_inc] = interpolate_vectors(x, y, ix, iy, L, 1);

hF = figure();
hA = axes(hF);
plot(x_ref_v_inc(1:80), 'LineWidth', 15, 'Color', [246 154 25]/255);
hold on;
plot(y_ref_v_inc(1:80), 'LineWidth', 15, 'Color', [25 172 246]/255);
hold off;
ylim([0 2.5]);
title(['Average Difference = ' num2str(dtw_d/length(ix))], 'FontSize', 15);
box off; grid off;
xticks([]);
set(get(hA, 'XAxis'), 'Visible', 'off');
set(get(hA, 'YAxis'), 'FontSize', 35, 'LineWidth', 5);

fig = gcf;
fig.Position = [100, 100, 1800, 700];

%% Helper Function: interpolate_vectors
function [x_ref_v_inc, y_ref_v_inc] = interpolate_vectors(x, y, ix, iy, L, gamma)
    % interpolate_vectors - Interpolates reference vectors based on DTW mapping.
    %
    % Syntax:
    %   [x_ref_v_inc, y_ref_v_inc] = interpolate_vectors(x, y, ix, iy, L, gamma)
    %
    % Inputs:
    %   x, y    - Original time course vectors (must be z-scored)
    %   ix, iy  - DTW alignment indices for x and y respectively
    %   L       - Length of the original time series (number of time points)
    %   gamma   - Exponent used to compute the reference difference (controls nonlinearity)
    %
    % Outputs:
    %   x_ref_v_inc, y_ref_v_inc - Interpolated reference vectors for x and y.
    %
    % The function computes an initial difference vector "v" based on the absolute
    % difference between the aligned values in x and y (with a sign based on the difference
    % in absolute amplitudes), raises it to the gamma power, and then uses PCHIP
    % interpolation to generate a continuous reference mapping.
    
    % Step 1: Compute the initial difference vector v (with sign)
    v = sign(abs(x(ix)) - abs(y(iy))) .* abs(x(ix) - y(iy)).^gamma;
    
    % Step 2: Define a sampling space matching the original number of time points
    x_samp = linspace(1, length(ix), L);
    
    % Step 3: Interpolate the difference vector using PCHIP interpolation
    v_int = interp1(v, x_samp, "pchip");

    % Step 4: Create reference vectors from the interpolated difference
    x_ref_v = v_int;
    y_ref_v = -v_int;
    
    % Enforce non-negativity of the reference vectors
    x_ref_v(x_ref_v < 0) = 0;
    y_ref_v(y_ref_v < 0) = 0;
    
    % Step 5: For this implementation, we return the computed vectors directly.
    x_ref_v_inc = x_ref_v;
    y_ref_v_inc = y_ref_v;
end
