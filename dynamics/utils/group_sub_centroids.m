function group_sub_clusters = group_sub_centroids(group_est, idx, cluster_num)
    % group_sub_centroids Groups sub-centroids into clusters based on provided indices
    %
    % This function calculates the cluster centroids for each subject and group
    % based on the given indices, grouping estimates by cluster number.
    %
    % INPUTS:
    %   group_est: A 3D numeric array of group estimates (subjects x time points x features).
    %   idx: A 2D numeric array of cluster indices (subjects x time points).
    %   cluster_num: An integer specifying the number of clusters.
    %
    % OUTPUT:
    %   group_sub_clusters: A 3D array of grouped sub-centroids (subjects x clusters x features).
    %
    % Example:
    %   group_sub_clusters = group_sub_centroids(group_est, idx, 5);

    % Input validation
    if nargin ~= 3
        error('Incorrect number of inputs. Expected 3 inputs: group_est, idx, and cluster_num.');
    end
    
    if ~isnumeric(group_est) || ~isnumeric(idx) || ~isnumeric(cluster_num)
        error('All inputs must be numeric.');
    end
    
    if ndims(group_est) ~= 3
        error('group_est must be a 3-dimensional array.');
    end
    
    if size(idx, 1) ~= size(group_est, 1)
        error('The number of rows in idx must match the number of subjects in group_est.');
    end
    
    % Initialize variables
    [sub_size, ~, feat_num] = size(group_est);
    group_sub_clusters = zeros(sub_size, cluster_num, feat_num); % Pre-allocate cluster centroids for subs for group1 

    % Computation loop
    for sub = 1:sub_size
        sub_group = idx(sub, :); % subject Group cluster index
        fprintf('sub: %d\n', sub);
        
        % Compute cluster centroid for each subject, each cluster, and each group        
        for kn = 1:cluster_num
            group_cluster_kn = find(sub_group == kn);

            if isscalar(group_cluster_kn)
                sub_cluster = group_est(sub, group_cluster_kn, :);
            elseif isempty(group_cluster_kn)
                sub_cluster = NaN([1, 1, feat_num]);
            else
                % Compute averaged FNCs for selected timepoints under cluster             
                sub_cluster = mean(group_est(sub, group_cluster_kn, :), 2);
            end
            
            group_sub_clusters(sub, kn, :) = squeeze(sub_cluster);
        end
    end
end