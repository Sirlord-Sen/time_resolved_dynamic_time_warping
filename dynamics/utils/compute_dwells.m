function dwell_count = compute_dwells(cluster_idx, kn_num)
% compute_dwells - Computes the number of dwells for a specified cluster.
%
% Syntax:
%   dwell_count = compute_dwells(cluster_idx, kn_num)
%
% Inputs:
%   cluster_idx : A vector containing cluster assignments over time.
%   kn_num      : The specific cluster number for which the dwell count is calculated.
%
% Outputs:
%   dwell_count : The number of dwells for the specified cluster.
%
% Description:
%   A "dwell" is defined as a contiguous sequence of time points where the cluster is active.
%   This function calculates the number of distinct dwell periods by identifying the start 
%   of each dwell (i.e., when the cluster becomes active after being inactive) in a vectorized manner.
%
% Example:
%   % Suppose cluster_idx is a vector of cluster assignments:
%   cluster_idx = [1 1 2 2 1 1 1 3 3 1];
%   % To compute the number of dwells for cluster 1:
%   dwell_count = compute_dwells(cluster_idx, 1);
%
% Author: Sir-Lord
% Date: [Today's Date]

    % Create a logical vector where the specified cluster is active.
    cluster_active = (cluster_idx == kn_num);
    
    % Identify the start of each dwell:
    % The first element is a dwell start if it is active.
    % A subsequent element is a dwell start if it is active and the previous element was inactive.
    dwell_starts = [cluster_active(1), cluster_active(2:end) & ~cluster_active(1:end-1)];
    
    % Sum the dwell start indicators to obtain the total dwell count.
    dwell_count = sum(dwell_starts);
end
