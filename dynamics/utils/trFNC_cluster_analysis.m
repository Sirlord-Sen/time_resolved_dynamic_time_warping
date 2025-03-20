function [mdt, fr, tm, tc] = trFNC_cluster_analysis(cluster_idx, cluster_num, exclude_self)
    % TRFNC_CLUSTER_ANALYSIS - Perform transition FNC (trFNC) analysis.
    % Computes various metrics such as mean dwell time (mdt), fraction of
    % recurrence (fr), and transition matrix (tm) for each subject.
    %
    % Inputs:
    %   cluster_idx : Matrix of cluster indices for each time point and subject.
    %                 Size: [num_subjects x num_timepoints]
    %   cluster_num : Number of clusters.
    %
    % Outputs:
    %   mdt : Mean dwell time for each subject and cluster. Size: [num_subjects x cluster_num]
    %   fr  : Fraction of recurrence for each subject and cluster. Size: [num_subjects x cluster_num]
    %   tm  : Transition probability matrix for each subject. Size: [num_subjects x cluster_num x cluster_num]
    %   tc  : Transition count matrix for each subject. Size: [num_subjects x cluster_num x cluster_num]
    
    if nargin < 3
        exclude_self = false; % Default to including self-transitions
    end

    % Check input arguments
    if nargin < 2
        error('Incorrect number of input arguments. Expected at least 2 inputs.');
    end
    
    if ~isnumeric(cluster_idx) || ~isnumeric(cluster_num)
        error('Inputs must be numeric.');
    end
    
    if size(cluster_idx, 2) < 2 || ~ismatrix(cluster_idx)
        error('cluster_idx must be a matrix with at least 2 columns.');
    end
    
    if ~isscalar(cluster_num) || cluster_num < 1 || mod(cluster_num, 1) ~= 0
        error('Invalid cluster_num. Must be a positive integer scalar.');
    end
    
    % Extract dimensions
    subs = size(cluster_idx, 1);
    tp = size(cluster_idx, 2);
    
    % Initialize output variables
    mdt = zeros(subs, cluster_num); %Mean-dwell time
    fr = zeros(subs, cluster_num); %fraction rate
    tm = zeros(subs, cluster_num, cluster_num); %transition matrix
    tc = zeros(subs, cluster_num, cluster_num); %transition counts
    
    % Loop through subjects
    for sub = 1:subs
        sub_idx = cluster_idx(sub, :);
    
        % Compute transition matrix
        [tm(sub, :, :), tc(sub, :, :)] = compute_transition_matrix(sub_idx, cluster_num, exclude_self);
    
        % Compute mean dwell time and fraction of recurrence
        for kn = 1:cluster_num
            dwell_count = compute_dwells(sub_idx, kn);
            if dwell_count > 0
                mdt(sub, kn) = sum(sub_idx == kn) / dwell_count;
            else
                mdt(sub, kn) = 0;  % Handle clusters with no dwells
            end
            fr(sub, kn) = sum(sub_idx == kn) / tp;
        end
    end
end