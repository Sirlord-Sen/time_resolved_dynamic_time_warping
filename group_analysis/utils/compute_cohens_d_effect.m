function d_val = compute_cohens_d_effect(fnc_vec, group)
% compute_cohens_d_effect - Computes Cohen's d effect size for each FNC feature between two groups.
%
% Syntax:
%   d_val = compute_cohens_d_effect(fnc_vec, group)
%
% Inputs:
%   fnc_vec - A matrix (n x m) where n is the number of observations (subjects)
%             and m is the number of FNC features.
%   group   - A vector (n x 1) containing the group assignments for each subject.
%             This vector should have exactly two unique values. The group with
%             the lowest label is taken as group1 and the group with the higher
%             label is taken as group2.
%
% Outputs:
%   d_val   - A row vector (1 x m) containing Cohen's d effect sizes for each FNC feature.
%
% Description:
%   The function splits the input FNC features according to the group assignment 
%   (with the lowest label as group1 and the highest as group2) and then computes 
%   Cohen's d effect size for each feature. Cohen's d is computed as:
%
%         d = (mean1 - mean2) / pooled_std,
%
%   where pooled_std is the pooled standard deviation of the two groups.
%
% Example:
%   % Suppose fnc_vec is a 100x50 matrix and group is a 100x1 binary vector:
%   d_val = compute_cohens_d_effect(fnc_vec, group);
%

    % Ensure groups are sorted in ascending order; group1 corresponds to the lowest label.
    groupings = sort(unique(group), "ascend");
    if numel(groupings) ~= 2
        error('The "group" input must contain exactly two unique values.');
    end
    
    % Assign group1 (lower label) and group2 (higher label)
    group1_idx = group == groupings(1);
    group2_idx = group == groupings(2);
    
    % Extract FNC features for each group.
    group1_fnc = fnc_vec(group1_idx, :);
    group2_fnc = fnc_vec(group2_idx, :);
    
    % Determine number of features.
    num_features = size(group1_fnc, 2);
    d_val = zeros(1, num_features);
    
    % Loop over each feature and compute Cohen's d.
    for k = 1:num_features
        x = group1_fnc(:, k);
        y = group2_fnc(:, k);
        
        % Compute group means.
        mean1 = mean(x);
        mean2 = mean(y);
    
        % Compute group standard deviations.
        std1 = std(x);
        std2 = std(y);
    
        % Get sample sizes.
        n1 = length(x);
        n2 = length(y);
    
        % Compute the pooled standard deviation.
        pooled_std = sqrt(((n1 - 1) * std1^2 + (n2 - 1) * std2^2) / (n1 + n2 - 2));
    
        % Compute Cohen's d.
        d_val(k) = (mean1 - mean2) / pooled_std;
    end
end
