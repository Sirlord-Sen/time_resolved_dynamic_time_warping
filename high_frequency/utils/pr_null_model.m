function surr_data = pr_null_model(data, num_surrogate)
% PR_NULL_MODEL Generate phase-randomized surrogate time series.
%
%   surr_data = PR_NULL_MODEL(data, num_surrogate) generates surrogate 
%   data by performing phase randomization on the input data. The 
%   function preserves the power spectrum of the original data while 
%   randomizing its phase. The number of surrogates to generate is 
%   specified by num_surrogate.
%
%   INPUTS:
%       data          - Input data time courses. A 3D array with dimensions 
%                       [subjects x time points x components].
%       num_surrogate - The number of surrogate data sets to generate.
%
%   OUTPUTS:
%       surr_data - Phase-randomized surrogate time series. A 4D array with 
%                   dimensions [subjects x time points x components x num_surrogate].
%
%   Example:
%       % Generate surrogate data for 10 subjects, 1000 time points, 
%       % 3 components, and 50 surrogates.
%       data = randn(10, 1000, 3); % Example input data
%       num_surrogate = 50;
%       surr_data = pr_null_model(data, num_surrogate);
%
%   Notes:
%       - This function handles both even and odd numbers of time points.
%       - The DC component (0 Hz) and, in the case of even time points, the 
%         Nyquist frequency are kept unchanged.
%
%   Author: [Sir-Lord]
%   Version: 1.0

    % Validate inputs
    if nargin < 2
        error('Not enough input arguments. Please provide both data and num_surrogate.');
    end

    if ~isnumeric(data) || ndims(data) ~= 3
        error('Input data must be a 3D numeric array of dimensions [subjects x time points x components].');
    end

    if ~isnumeric(num_surrogate) || ~isscalar(num_surrogate) || num_surrogate <= 0
        error('Input num_surrogate must be a positive scalar integer.');
    end

    % Get dimensions
    [subs, tps, comps] = size(data);
    
    % Compute half-length of the time points
    len_h = floor(tps/2);
    
    % Define frequency intervals based on whether tps is even or odd
    if mod(tps, 2) == 0
        % Even number of time points
        intv_1 = 2:len_h;
    else
        % Odd number of time points
        intv_1 = 2:len_h+1;
    end
    intv_2 = len_h + 2 : tps;
    
    % Initialize output array for surrogate time series
    surr_data = zeros(subs, tps, comps, num_surrogate);

    % Generate surrogate data
    for k = 1:num_surrogate
        for sub = 1:subs
            % Generate random phases for shuffling
            ph_rand = rand(length(intv_1), 1);
            
            % Fourier transform of the original data for this subject
            data_fft = fft(squeeze(data(sub, :, :))); 

            % Create random phases for positive frequencies
            ph_intv_1 = repmat(exp(2*pi*1i*ph_rand), 1, comps);
    
            % Apply conjugate symmetry for negative frequencies
            ph_intv_2 = conj(flipud(ph_intv_1));
    
            % Apply phase randomization
            surr_fft = data_fft;
            surr_fft(intv_1, :) = data_fft(intv_1, :) .* ph_intv_1;
            surr_fft(intv_2, :) = data_fft(intv_2, :) .* ph_intv_2;
    
            % Inverse Fourier transform to get the surrogate time series
            surr_data(sub, :, :, k) = real(ifft(surr_fft));
        end
    end

end
