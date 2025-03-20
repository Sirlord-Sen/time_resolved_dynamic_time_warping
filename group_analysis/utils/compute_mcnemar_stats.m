function [p, s] = compute_mcnemar_stats(x, y)
% compute_mcnemar_stats - Computes McNemar's test statistics for two binary vectors.
%
% Syntax:
%   [p, s] = compute_mcnemar_stats(x, y)
%
% Inputs:
%   x - A binary vector (e.g., logical or numeric vector with values 0 and 1) 
%       representing significance status for condition x.
%   y - A binary vector representing significance status for condition y.
%
% Outputs:
%   p - The p-value computed from McNemar's test based on the 2x2 contingency table.
%   s - The sign of the difference in the total number of significant features 
%       (i.e., sign(sum(x) - sum(y))).
%
% Description:
%   This function creates a 2x2 contingency table using the input vectors x and y,
%   where:
%       a = number of features significant in both x and y,
%       b = number of features significant in x but not in y,
%       c = number of features significant in y but not in x,
%       d = set equal to a (symmetry assumption).
%
%   It then computes the p-value using McNemar's test (via the mcnemar function)
%   and determines the direction of difference (s) between the sums of x and y.
%
% Example:
%   x = [1 0 1 1 0];
%   y = [1 1 0 1 0];
%   [p, s] = compute_mcnemar_stats(x, y);
%

    % Create a 2x2 contingency table.
    s_data = zeros(2, 2);
    common_count = numel(intersect(find(x), find(y)));
    s_data(1,1) = common_count;
    s_data(2,2) = common_count;
    s_data(1,2) = numel(setdiff(find(x), find(y)));
    s_data(2,1) = numel(setdiff(find(y), find(x)));
    
    % Compute the p-value using McNemar's test.
    p = mcnemar(s_data);
    
    % Record the direction of difference (sign of the difference in the number of significant features).
    s = sign(sum(x) - sum(y));
end
