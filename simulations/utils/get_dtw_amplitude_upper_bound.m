function [dtw_max_dist, dtw_norm_max_dist] = get_dtw_amplitude_upper_bound(x, win_size, gamma, F, Fs, amp_max, num_trials)
% get_dtw_amplitude_upper_bound - Calculate the upper bound of DTW distances between a sine wave 
% and a randomly amplitude-perturbed version using multiple trials.
%
% Syntax:
%   [dtw_max_dist, dtw_norm_max_dist] = get_dtw_amplitude_upper_bound(x, win_size, gamma, F, Fs, amp_max, num_trials)
%
% Parameters:
%   x           - Original sine wave (vector).
%   win_size    - DTW warping window size constraint (scalar). Determines the range of indices 
%                 considered during DTW computation.
%   gamma       - Exponent for the distance measure in DTW (scalar, default is 1 if not provided).
%   F           - Frequency of the sine wave (Hz). (Note: F is not used directly in this function,
%                 but is included for consistency with other functions.)
%   Fs          - Sampling frequency (Hz).
%   amp_max     - Maximum amplitude deviation for the random amplitude perturbation.
%   num_trials  - Number of trials to perform for the random amplitude perturbation (integer).
%
% Returns:
%   dtw_max_dist       - Array of DTW distances for each trial.
%   dtw_norm_max_dist  - Array of normalized DTW distances for each trial (DTW distance divided by path length).

    % Time vector based on the sampling frequency and length of x
    Ts = 1/Fs;
    L = length(x);
    t = (0:Ts:(Ts*L - Ts))'; 
    
    % Pre-allocate arrays for DTW distances
    dtw_max_dist = zeros(1, num_trials);
    dtw_norm_max_dist = zeros(1, num_trials);

    % Loop to calculate DTW distances for each random amplitude perturbation trial
    for i = 1:num_trials
        % Create random amplitude noise between 0 and amp_max
        random_amp_noise = amp_max * rand(size(t));
        
        % Generate the second sine wave with random amplitude perturbation
        y = random_amp_noise .* x;
        
        % Calculate DTW distance using the dtw_custom function
        [~, dtw_max_dist(i), ix] = dtw_custom(zscore(x), zscore(y), win_size, gamma);
        
        % Normalize DTW distance by the length of the warping path
        dtw_norm_max_dist(i) = dtw_max_dist(i) / length(ix);
    end
end