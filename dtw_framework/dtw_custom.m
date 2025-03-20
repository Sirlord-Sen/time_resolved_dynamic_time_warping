function [D, cost, ix, iy] = dtw_custom(x, y, win_size, gamma)
% dtw_custom - Compute the Dynamic Time Warping (DTW) distance and optimal warping path
% between two time series.
%
% Syntax:
%   [D, cost, ix, iy] = dtw_custom(x, y)
%   [D, cost, ix, iy] = dtw_custom(x, y, win_size)
%   [D, cost, ix, iy] = dtw_custom(x, y, win_size, gamma)
%
% Inputs:
%   x        - First time series (vector).
%   y        - Second time series (vector).
%   win_size - (Optional) Warping window constraint. If not provided or empty,
%              it is set to the length of the shorter time series (or 1 if both
%              signals are of equal length).
%   gamma    - (Optional) Exponent for the distance measure. Default is 1.
%
% Outputs:
%   D        - Accumulated cost matrix of size (n+1)x(m+1), where n=length(x) and m=length(y).
%   cost     - Final DTW cost (scalar), located at D(n+1, m+1).
%   ix       - Indices in x corresponding to the optimal warping path.
%   iy       - Indices in y corresponding to the optimal warping path.
%
% Example:
%   x = randn(100,1);
%   y = randn(100,1);
%   [D, cost, ix, iy] = dtw_custom(x, y); % Uses default win_size and gamma
%
% Author: Sir-Lord Wiafe

    % Determine lengths of the input time series
    n = length(x);
    m = length(y);

    % Set default window size if not provided or empty
    if nargin < 3 || isempty(win_size)
        if n == m
            win_size = m;
        else
            win_size = min(n, m);
        end
    end

    % Set default gamma if not provided or empty
    if nargin < 4 || isempty(gamma)
        gamma = 1;
    end

    % Ensure the window size is at least the absolute difference in lengths
    win_size = max(win_size, abs(n - m));

    % Initialize the accumulated cost matrix with infinities.
    % The matrix has dimensions (n+1)x(m+1) to accommodate boundary conditions.
    D = inf(n+1, m+1);
    D(1, 1) = 0;
    
    % Define the distance function using the provided gamma exponent.
    % For gamma=1, this is equivalent to the absolute difference.
    distFun = @(a, b) abs(a - b)^gamma;
    
    % Populate the cost matrix with the dynamic programming approach.
    for i = 1:n
        % Define valid index range for y based on win_size constraint.
        jStart = max(1, i - win_size);
        jEnd = min(m, i + win_size);
        for j = jStart:jEnd
            cost_val = distFun(x(i), y(j));
            D(i+1, j+1) = cost_val + min([D(i, j+1), D(i+1, j), D(i, j)]);
        end
    end
    
    % The final DTW cost is found in the bottom-right corner of the matrix.
    cost = D(n+1, m+1);
    
    % Backtrack to obtain the optimal warping path.
    [ix, iy] = traceback(D, n, m);
end

function [ix_out, iy_out] = traceback(D, n, m)
% traceback - Backtrack through the accumulated cost matrix to extract the
% optimal warping path for DTW.
%
% Syntax:
%   [ix_out, iy_out] = traceback(D, n, m)
%
% Inputs:
%   D - Accumulated cost matrix from dtw_custom (size: (n+1)x(m+1)).
%   n - Length of the first time series.
%   m - Length of the second time series.
%
% Outputs:
%   ix_out - Indices in the first time series corresponding to the optimal path.
%   iy_out - Indices in the second time series corresponding to the optimal path.
%
% Example:
%   [ix, iy] = traceback(D, n, m);

    % Start at the bottom-right corner of the cost matrix.
    i = n + 1;
    j = m + 1;
    
    % Initialize arrays to hold the warping path indices.
    ix = [];
    iy = [];
    
    % Backtrack until reaching the first row or column.
    while i > 1 && j > 1
        ix = [i - 1; ix];
        iy = [j - 1; iy];
        
        % Handle edge cases when near the matrix boundaries.
        if i == 2
            j = j - 1;
        elseif j == 2
            i = i - 1;
        else
            % Choose the direction that minimizes the accumulated cost.
            [~, min_idx] = min([D(i-1, j-1), D(i-1, j), D(i, j-1)]);
            if min_idx == 1
                i = i - 1;
                j = j - 1;
            elseif min_idx == 2
                i = i - 1;
            else
                j = j - 1;
            end
        end
    end
    
    % If the path hasn't reached the first column, add the remaining steps.
    while i > 1
        ix = [i - 1; ix];
        iy = [1; iy];
        i = i - 1;
    end
    
    % If the path hasn't reached the first row, add the remaining steps.
    while j > 1
        ix = [1; ix];
        iy = [j - 1; iy];
        j = j - 1;
    end
    
    % Output the optimal warping path indices.
    ix_out = ix(2:end);
    iy_out = iy(2:end);
end
