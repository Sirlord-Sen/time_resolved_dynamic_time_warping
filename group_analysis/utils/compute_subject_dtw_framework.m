function [full_dtw_d, full_dtw_d_norm, full_tr_dtw, full_d_dtw] = compute_subject_dtw_framework(subTcs, gamma_vec, Tr, band)
% compute_subject_dtw_framework - Computes DTW metrics for all subjects and unique component pairs.
%
% Syntax:
%   [full_dtw_d, full_dtw_d_norm, full_tr_dtw, full_d_dtw] = compute_group_dtw(subTcs, gamma_vec, Tr, band)
%
% Inputs:
%   subTcs    - A 3D matrix of time courses with dimensions 
%               [num_subjects x num_timepoints x num_components].
%   gamma_vec - A vector of gamma values to use in the DTW computation.
%   Tr        - Repetition time (sampling interval) in seconds.
%   band      - A two-element vector specifying the frequency band for DTWFramework.
%
% Outputs:
%   full_dtw_d      - A 3D array of raw DTW distances with dimensions 
%                     [num_gammas x num_subjects x num_features], where num_features = nchoosek(num_components, 2).
%   full_dtw_d_norm - A 3D array of normalized DTW distances with the same dimensions.
%   full_tr_dtw     - A 4D array of time-resolved DTW distances with dimensions 
%                     [num_gammas x num_subjects x num_timepoints x num_features].
%   full_tr_dtw     - A 4D array of directional DTW distances with same dimensions as full_tr_dtw 
%
% Description:
%   This function computes DTW distances between all unique pairs of components for each subject 
%   and for each gamma value provided. For each subject, it computes a lower-triangular matrix of 
%   DTW distances between the components, vectorizes the unique values using icatb_mat2vec, and 
%   stores the result. The DTWFramework function is used to compute both the raw and normalized 
%   DTW distances.
%
% Example:
%   % Assume subTcs is a 3D matrix of dimensions [100 subjects x 150 timepoints x 10 components],
%   % gamma_vec = 0.25:0.25:3, Tr = 0.72, and band = [0.01 0.198]:
%   [full_dtw, full_dtw_norm] = compute_group_dtw(subTcs, 0.25:0.25:3, 0.72, [0.01 0.198]);

    % Determine dimensions from input
    num_subjects   = size(subTcs, 1);
    num_timepoints = size(subTcs, 2);
    num_components = size(subTcs, 3);
    num_gammas     = length(gamma_vec);
    num_features   = nchoosek(num_components, 2);
    
    % Preallocate output matrices
    full_dtw_d      = zeros(num_gammas, num_subjects, num_features);
    full_dtw_d_norm = zeros(num_gammas, num_subjects, num_features);
    full_tr_dtw = zeros(num_gammas, num_subjects, num_timepoints, num_features);
    full_d_dtw = zeros(num_gammas, num_subjects, num_timepoints, num_features);
    
    % Loop over each gamma value
    for gamma_idx = 1:num_gammas
        gamma = gamma_vec(gamma_idx);
        
        for sub = 1:num_subjects
            % Display progress information in the command window.
            fprintf('Processing gamma = %.2f (%d/%d), subject = %d/%d\n', ...
                gamma, gamma_idx, num_gammas, sub, num_subjects);
            
            % Initialize matrices to store DTW distances for the current subject.
            % These matrices are of size [num_components x num_components] where
            % only the lower triangular part (unique component pairs) is filled.
            sub_dtw_d      = zeros(num_components, num_components);
            sub_dtw_d_norm = zeros(num_components, num_components);
            sub_tr_dtw     = zeros(num_timepoints, num_components, num_components);
            sub_d_dtw      = zeros(num_timepoints, num_components, num_components);
            
            % Loop over all unique component pairs (comp1 > comp2)
            for comp1 = 2:num_components
                for comp2 = 1:(comp1 - 1)
                    % Extract the time courses for the two components and z-score them.
                    x = squeeze(subTcs(sub, :, comp1));
                    y = squeeze(subTcs(sub, :, comp2));
                    x = zscore(x);
                    y = zscore(y);
                    
                    % Compute DTW metrics using the custom DTW framework.
                    [dtw_d, dtw_d_norm, tr_dtw, d_dtw] = DTWFramework(x, y, gamma, Tr, band);
                    
                    % Store the computed DTW distances in the lower-triangular matrix.
                    sub_dtw_d(comp1, comp2) = dtw_d;
                    sub_dtw_d_norm(comp1, comp2) = dtw_d_norm;
                    sub_tr_dtw(:, comp1, comp2) = tr_dtw;
                    sub_d_dtw(:, comp1, comp2) = d_dtw;
                end
            end
            
            % Vectorize the lower-triangular matrices and store them in the output arrays.
            full_dtw_d(gamma_idx, sub, :) = icatb_mat2vec(sub_dtw_d);
            full_dtw_d_norm(gamma_idx, sub, :) = icatb_mat2vec(sub_dtw_d_norm);
            full_tr_dtw(gamma_idx, sub, :, :) = icatb_mat2vec(sub_tr_dtw);
            full_d_dtw(gamma_idx, sub, :, :) = icatb_mat2vec(sub_d_dtw);
        end
    end
end
