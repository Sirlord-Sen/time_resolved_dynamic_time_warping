function win_size = calculate_window_size(TR, band)
    % calculate_window_size - Computes an odd window size for a given TR and frequency band.
    % 
    % Syntax:
    %   win_size = calculate_window_size(TR, band)
    %
    % Inputs:
    %   TR   - Sampling interval (in seconds). Must be a positive scalar.
    %   band - Frequency band specified as [minFreq maxFreq].
    %          Each frequency value must be positive.
    %
    % Output:
    %   win_size - Calculated window size, rounded up to the nearest odd integer.
    %
    % Example:
    %   win_size = calculate_window_size(2, [0.01 0.1]);

    % Input validation
    if nargin ~= 2
        error('calculate_window_size requires exactly two input arguments: TR and band.');
    end
    
    % Check TR is a positive scalar
    if ~isscalar(TR) || ~isnumeric(TR) || TR <= 0
        error('TR must be a positive scalar.');
    end

    % Check band is a vector of one or two positive values
    if ~isnumeric(band) || ~isvector(band) || isempty(band) || length(band) > 2 || any(band <= 0)
        error('band must be a vector containing one or two positive numeric values.');
    end
    
    % Use the first element of band as the primary frequency (minFreq or centerFreq)
    freq = band(1);

    % Calculate the initial window size using the formula
    win_size = sqrt(((0.88 / freq) * (1 / TR))^2 + 1);

    % Round up to the nearest integer and make it odd
    win_size = ceil(win_size); % Round up
    if mod(win_size, 2) == 0
        win_size = win_size + 1; % Ensure it's odd
    end
end
