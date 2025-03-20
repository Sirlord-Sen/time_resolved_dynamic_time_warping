function [transition_prob, transition_cnts] = compute_transition_matrix(sub_idx, kn_num, exclude_self)
    % COMPUTE_TRANSITION_MATRIX - Calculate transition probabilities for a given subject.
    % Inputs:
    %   sub_idx : Index vector representing cluster assignments over time for a subject.
    %   kn_num  : Number of clusters.
    %   exclude_self : Boolean flag to exclude self-transitions (default: false).
    % Outputs:
    %   transition_prob : Transition probabilities between clusters.
    %   transition_cnts : Transition counts between clusters.
    
    if nargin < 3
        exclude_self = false; % Default to including self-transitions
    end
    
    % Vectorized transition counting
    transitions = sub_idx(1:end-1) + (sub_idx(2:end) - 1) * kn_num;
    transition_cnts = accumarray(transitions', 1, [kn_num * kn_num, 1]);
    transition_cnts = reshape(transition_cnts, kn_num, kn_num);
    
    % Exclude self-transitions if requested
    if exclude_self
        transition_cnts(logical(eye(kn_num))) = 0; % Set diagonal elements to zero
    end
    
    % Normalize to compute transition probabilities
    transition_prob = transition_cnts ./ sum(transition_cnts, 2);
    transition_prob(isnan(transition_prob)) = 0;  % Handle divisions by zero
end
