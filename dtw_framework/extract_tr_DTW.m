function [tr_dtw, d_dtw_x, d_dtw_y] = extract_tr_DTW(x, y, ix, iy, L, gamma)
% extract_tr_DTW - Extracts the time-resolved amplitude mismatch between two signals after DTW warping.
%
% Syntax:
%   [tr_dtw, d_dtw_x, d_dtw_y] = extract_tr_DTW(x, y, ix, iy, L, gamma)
%
% Inputs:
%   x, y    - Original time course vectors (should be z-scored).
%   ix, iy  - DTW alignment indices for x and y, respectively.
%   L       - Desired length for the interpolated vectors (number of time points).
%   gamma   - Exponent to control nonlinearity when computing the difference vector.
%
% Outputs:
%   tr_dtw   - Absolute value of the interpolated difference vector, representing the
%              time-resolved amplitude mismatch between the two signals.
%   d_dtw_x  - Interpolated tr_DTW with x as reference (non-negative).
%   d_dtw_y  - Interpolated tr_DTW with y as reference (non-negative).
%
% Example:
%   [tr, d_x, d_y] = extract_tr_DTW(x, y, ix, iy, 1000, 1);
%

    % Step 1: Compute the initial difference vector v with sign.
    v = sign(abs(x(ix)) - abs(y(iy))) .* abs(x(ix) - y(iy)).^gamma;
    
    % Step 2: Define the sampling space matching the desired length L.
    x_samp = linspace(1, length(ix), L);
    
    % Step 3: Interpolate the difference vector v using PCHIP interpolation.
    v_int = interp1(v, x_samp, 'pchip');
    
    % Step 4: Create reference vectors for x and y.
    d_dtw_x = v_int;
    d_dtw_y = -v_int;
    
    % Enforce non-negativity of the reference vectors.
    d_dtw_x(d_dtw_x < 0) = 0;
    d_dtw_y(d_dtw_y < 0) = 0;
    
    % Step 5: The absolute interpolated difference represents the time-resolved mismatch.
    tr_dtw = abs(v_int);
end
