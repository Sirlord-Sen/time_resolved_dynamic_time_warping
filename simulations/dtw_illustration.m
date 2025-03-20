%% DTW Simulation and Visualization
% This script simulates time warping on a signal and demonstrates the alignment
% using Dynamic Time Warping (DTW). Two signals are generated:
%   1. A base signal created by bandpass filtering and normalization.
%   2. A time-warped version of the base signal via interpolation.
%
% The DTW alignment is then computed and visualized.
%
%
% Author: [Sir-Lord Wiafe]
% Date: [20th February 2025]
% -------------------------------------------------------------------------

%% Clear Environment and Set Paths
clc; clear;
addpath(genpath('./utils/'))
%% Simulation Parameters
Tr = 1;                          % Sampling interval (sec)
Fs = 1/Tr;                       % Sampling frequency (Hz)
L = 200;                         % Signal length (number of samples)
timeVector = (0:1/Fs:(Tr*L - 1/Fs))';  % Time vector

% Warping parameters
max_stretch = 15;                
warp_factor = max_stretch / Tr;  % Amplitude of time modulation

%% Generate Time Warping Function
% Create a time warping function using bandpass-filtered noise.
cutoff = [0.5 1.5];
band = [0.01 0.02];
warp_func = bandpass_filtering(randn(L, 1), Tr, band, cutoff);
warp_func = warp_func .* (warp_factor / max(abs(warp_func)));
time_warp = timeVector + warp_func;

figure;
plot(warp_func, '--k', LineWidth=3)
title('Warping function');
xlabel('Time Index');
ylabel('Amplitude');
grid on;


%% Generate Base Signal and Create Warped Signal
band = [0.01 0.1];
cutoff = [0.5 1.5];
% Generate the base signal using bandpass filtering and normalization
signal1 = bandpass_filtering(randn(L, 1), Tr, band, cutoff);
signal1 = zscore(signal1);
signal1 = (signal1 / max(abs(signal1))) * 0.2;

% Create the warped signal via interpolation
signal2 = interp1(timeVector, signal1, time_warp, 'linear', 'extrap');
signal2 = zscore(signal2);
signal1 = zscore(signal1);

%% Plot Original and Warped Signals
figure;
plot(signal1, 'LineWidth', 5, 'Color', [98 75 131]/255);
hold on;
plot(signal2, 'LineWidth', 5, 'Color', [158 151 36]/255);
hold off;
title('Original and Warped Signals');
xlabel('Time Index');
ylabel('Normalized Amplitude');
legend('Original signal', 'Warped signal');
grid on;

%% Compute DTW Alignment
[~, ix, iy] = dtw(signal1, signal2);

%% Plot Detailed DTW Alignment with Vertical Shift
plotDTWAlignment(signal1, signal2, ix, iy, 7);
set(gcf, 'Position', [100, 100, 1200, 700]);

%% Function Definitions
function plotDTWAlignment(x, y, ix, iy, shift)
% plotDTWAlignment - Visualizes DTW alignment between two time series.
%
% Syntax:
%   plotDTWAlignment(x, y, ix, iy, shift)
%
% Inputs:
%   x     - First time series (vector).
%   y     - Second time series (vector).
%   ix    - Indices from the first time series mapping (vector).
%   iy    - Indices from the second time series mapping (vector).
%   shift - Vertical offset applied to the second time series.
%
% Example:
%   plotDTWAlignment(signal1, signal2, ix, iy, 7);
%
% Author: [Your Name]
% Date: [Today's Date]

    % Validate input arguments
    if nargin ~= 5
        error('plotDTWAlignment:InvalidInput', 'Function requires exactly 5 inputs.');
    end

    if ~isvector(x) || ~isvector(y)
        error('plotDTWAlignment:InvalidInput', 'Inputs x and y must be vectors.');
    end

    if length(ix) ~= length(iy)
        error('plotDTWAlignment:InvalidInput', 'Mapping indices ix and iy must have the same length.');
    end

    if ~isscalar(shift)
        error('plotDTWAlignment:InvalidInput', 'Input shift must be a scalar.');
    end

    % Create a new figure for alignment visualization
    figure;
    hold on;
    
    % Plot the first time series
    plot(x, 'LineWidth', 7, 'Color', [98 75 131]/255);
    
    % Plot the second time series with a vertical shift for clarity
    plot(y - shift, 'LineWidth', 7, 'Color', [158 151 36]/255);
    
    % Draw alignment lines between the matched indices
    for idx = 1:length(ix)
        plot([ix(idx) iy(idx)], [x(ix(idx)) (y(iy(idx)) - shift)], 'k', 'LineWidth', 1);
    end
    
    % Remove axis labels and ticks for a cleaner visualization
    xticks([]);
    yticks([]);
    ax = gca;
    ax.XColor = 'none';
    ax.YColor = 'none';
    hold off;
end

