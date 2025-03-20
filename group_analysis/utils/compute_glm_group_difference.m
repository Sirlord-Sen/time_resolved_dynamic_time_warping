function [q_values, p_values, t_values] = compute_glm_group_difference(fnc_vec, group, covariates)
% compute_glm_group_difference - Computes group differences across FNC features using GLM.
%
% Syntax:
%   [q_values, p_values, t_values] = compute_glm_group_difference(fnc_matrix, group, covariates)
%
% Inputs:
%   fnc_vec - A matrix (n x m), where n is the number of observations (e.g., subjects)
%                and m is the number of FNC features. Each column represents an FNC feature.
%
%   group      - A binary vector (n x 1) representing the group assignment of each subject.
%                Typically, values are 0 and 1, where 1 might represent a specific group (e.g., patients)
%                and 0 another (e.g., controls).
%
%   covariates - A table (n x p) of covariates for each subject, where n is the number of
%                observations and p is the number of covariates. This table can include any
%                relevant variables, such as age, sex, site, or other covariates.
%
% Outputs:
%   q_values   - A vector (m x 1) of q-values for each FNC feature after applying
%                Benjamini-Hochberg FDR correction.
%
%   p_values   - A vector (m x 1) of p-values for each FNC feature before FDR correction.
%
%   t_values   - A vector (m x 1) of t-statistics for the group effect for each FNC feature.
%
% Example:
%   % Example Data (Randomized for demonstration)
%   n = 311;          % Number of subjects
%   m = 100;          % Number of FNC features
%   fnc_vec = randn(n, m);       % Example FNC matrix (311 subjects, 100 features)
%   group = randi([0, 1], n, 1);      % Random binary group assignment
%   
%   % Creating a table of covariates
%   age = randi([20, 80], n, 1);                % Random ages between 20 and 80
%   sex = randi([0, 1], n, 1);                  % Random binary sex assignment
%   site = randi([1, 3], n, 1);                 % Random site assignment with 3 possible sites
%   mfd = randn(n, 1);                          % Random mean frame displacement assignment
%   covariates = table(age, sex, site, mfd);    % Table of covariates
%   
%   % Call the function
%   [q_values, p_values, t_values] = compute_glm_group_difference(fnc_matrix, group, covariates);
%

    % Validate that the number of observations in covariates matches fnc_matrix
    if height(covariates) ~= size(fnc_vec, 1)
        error('The number of rows in covariates must match the number of rows in fnc_matrix.');
    end
    
    % Determine the number of FNC features
    num_features = size(fnc_vec, 2);
    
    % Preallocate vectors for p-values and t-values
    p_values = zeros(num_features, 1);
    t_values = zeros(num_features, 1);
    
    % Add 'group' as the first variable in the covariates table
    covariates = addvars(covariates, group, 'NewVariableNames', 'group', 'Before', 1);
    
    % Loop through each FNC feature and fit a GLM
    for feature_idx = 1:num_features
        % Extract the current FNC feature data (a column vector)
        fnc_feature = fnc_vec(:, feature_idx);
        
        % Add the FNC feature as the response variable to the covariates table
        data = addvars(covariates, fnc_feature, 'NewVariableNames', 'fnc');
        
        % Fit a GLM with 'group' as the predictor and other covariates as controls.
        mdl = fitglm(data, 'ResponseVar', 'fnc');
        
        % Extract the p-value and t-statistic for the 'group' effect.
        % (The row with 'group' is indexed by its name in the table.)
        p_values(feature_idx) = mdl.Coefficients{'group','pValue'};
        t_values(feature_idx) = mdl.Coefficients{'group','tStat'};
    end
    % Apply Benjamini-Hochberg FDR correction to obtain q-values.
    q_values = mafdr(p_values, 'BHFDR', true);
end
