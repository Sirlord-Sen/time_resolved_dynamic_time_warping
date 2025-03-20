function filteredSignal = bandpass_filtering(signal, Tr, band, cutoff)
% bandpass_filtering - Applies a Butterworth bandpass filter to a time series.
%
% Syntax:
%   filteredSignal = bandpass_filtering(Tr, band, signal)
%
% Inputs:
%   signal - Input time series (vector).
%   Tr     - Sampling interval (seconds).
%   band   - Two-element vector specifying the [low high] passband frequencies (Hz).
%   cutoff - Two-element vector specifying the factor of passband.

%
% Optional Inputs:
%   pad_size - Number of padding points to mitigate edge effects (default: 100).
%
% Output:
%   filteredSignal - The bandpass-filtered time series.
%
% Example:
%   filteredSignal = bandpass_filtering(1, [0.01 0.1], randn(200,1));
%
% Author: Sir-Lord

    % Default padding size (if required for further processing)
    if nargin < 5
        pad_size = 100;
    end
    
    % Define sampling and Nyquist frequencies
    Fs = 1 / Tr;
    nyquist = Fs / 2;
    
    % Normalize the passband and define a wider stopband for robust filtering
    Wp = band / nyquist;
    Ws = [band(1)*cutoff(1) band(2)*cutoff(2)]/nyquist;
    
    % Define filter design parameters
    Rp = 3;   % Passband ripple (dB)
    Rs = 30;  % Stopband attenuation (dB)
    
    % Design a Butterworth bandpass filter
    [n, Wn] = buttord(Wp, Ws, Rp, Rs);
    [b, a] = butter(n, Wn, 'bandpass');
    
    % Apply zero-phase filtering to avoid phase distortion
    filteredSignal = filtfilt(b, a, signal);
end
