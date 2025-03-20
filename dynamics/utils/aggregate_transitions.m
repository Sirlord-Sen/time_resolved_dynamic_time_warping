function aggregate_transition_matrix = aggregate_transitions(transition_cnt)
% aggregate_transitions - Aggregates transition count matrices and computes an 
%                           aggregate transition probability matrix.
%
% Syntax:
%   aggregate_transition_matrix = aggregate_transitions(transition_cnt)
%
% Inputs:
%   transition_cnt - A 3D array (subjects x num_states x num_states) where each
%                    slice corresponds to a transition count matrix for a subject.
%
% Outputs:
%   aggregate_transition_matrix - A matrix (num_states x num_states) representing the
%                                 aggregate transition probability matrix. This matrix
%                                 is computed by summing transition counts across subjects
%                                 and normalizing each row. If a row sums to zero, uniform
%                                 probabilities are assigned.
%
% Description:
%   This function aggregates the transition count matrices from multiple subjects by
%   summing them and then normalizes each row to obtain a transition probability matrix.
%   In cases where a row sums to zero, the function assigns a uniform distribution across
%   all states and issues a warning.
%
% Author: [Sir-Lord]
%
% Example:
%   % Given transition_cnt as a 3D array of transition counts:
%   agg_matrix = aggregate_transitions(transition_cnt);

    % Determine the number of states from the transition count matrix.
    num_states = size(transition_cnt, 2);
    
    % Initialize aggregate count matrix.
    aggregate_counts = zeros(num_states, num_states);
    
    % Sum transitions across all subjects.
    for subj = 1:size(transition_cnt, 1)
        aggregate_counts = aggregate_counts + squeeze(transition_cnt(subj, :, :));
    end
    
    % Initialize aggregate transition probability matrix.
    aggregate_transition_matrix = zeros(num_states, num_states);
    
    % Normalize each row of the aggregate counts to create a probability matrix.
    for from = 1:num_states
        row_sum = sum(aggregate_counts(from, :));
        if row_sum > 0
            aggregate_transition_matrix(from, :) = aggregate_counts(from, :) / row_sum;
        else
            % Handle zero rows by assigning uniform probabilities.
            aggregate_transition_matrix(from, :) = ones(1, num_states) / num_states;
            warning(['Row ', num2str(from), ' sums to zero. Assigned uniform probabilities.']);
        end
    end
end
