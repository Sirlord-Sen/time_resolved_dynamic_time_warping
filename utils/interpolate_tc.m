function int_tc = interpolate_tc(tc, desired_len, method)
% interpolate_tc - Resample a time course vector to a desired length.
%
% Syntax:
%   int_tc = interpolate_tc(tc, desired_len, method)
%
% Inputs:
%   tc          - Input time course vector.
%   desired_len - The desired length (number of samples) after interpolation.
%   method      - Interpolation method (e.g., 'linear', 'pchip', etc.).
%
% Outputs:
%   int_tc      - The interpolated time course vector.
%
% Example:
%   % Resample a time course vector 'tc' to 500 samples using pchip interpolation.
%   new_tc = interpolate_tc(tc, 500, 'pchip');
%
% Author: Sir-Lord

    tc_len = length(tc);
    % Create a vector of new sample indices linearly spaced from 1 to tc_len
    new_index = linspace(1, tc_len, desired_len);
    % Interpolate using the original indices (1:tc_len) and the specified method
    int_tc = interp1(1:tc_len, tc, new_index, method);
end
