function [q_values, p_values, t_values, b_values] = compute_glm_associations(fnc, covariates, clinical_score, dist)
% compute_glm_associations - Computes GLM statistics for each FNC feature using a combined covariates matrix.
%
% Syntax:
%   [q_values, p_values, t_values, b_values] = compute_glm_feature_stats(fnc, covariates, clinical_score, dist)
%
% Inputs:
%   fnc            - A matrix (n x m) of FNC features, where n is the number of subjects 
%                    and m is the number of features.
%   covariates     - A matrix (n x p) of covariate predictors for each subject. This should include
%                    all predictors that are constant across features (e.g., mfd, age, sex, and the
%                    site dummies, with the first column of the site dummies already excluded).
%   clinical_score - A vector (n x 1) representing the clinical scores (response variable).
%   dist           - (Optional) Distribution for the GLM. Default is 'normal'. Supported options include
%                    'normal' and 'poisson'.
%
% Outputs:
%   q_values  - A row vector (1 x m) of q-values for each FNC feature after FDR correction.
%   p_values  - A row vector (1 x m) of uncorrected p-values for the FNC feature predictor.
%   t_values  - A row vector (1 x m) of t-statistics for the FNC feature predictor.
%   b_values  - A row vector (1 x m) of beta coefficients for the FNC feature predictor.
%
% Description:
%   The function fits a GLM for each FNC feature by using the FNC feature as the first predictor 
%   (across subjects) and combining it with the provided covariates. The effect of the FNC feature, 
%   corresponding to the second coefficient (after the intercept), is extracted. The function returns 
%   the p-value, t-statistic, and beta coefficient for the FNC feature, and applies Benjamini-Hochberg 
%   FDR correction to the p-values.
%
% Example:
%   % Suppose fnc is a 100x50 matrix, covariates is a 100x4 matrix (e.g., [mfd, age, sex, site_dummy]),
%   % and clinical_score is a 100x1 vector:
%   [q_vals, p_vals, t_vals, b_vals] = compute_glm_feature_stats(fnc, covariates, clinical_score, 'normal');

    % Set default distribution to 'normal' if not provided.
    if nargin < 4 || isempty(dist)
        dist = 'normal';
    end

    % Determine the link function based on the distribution.
    if strcmp(dist, 'poisson')
        link = 'log';
    elseif strcmp(dist, 'normal')
        link = 'identity';
    else
        error('Unsupported distribution: %s. Supported distributions are ''normal'' and ''poisson''.', dist);
    end

    % Determine the number of features.
    num_features = size(fnc, 2);

    % Preallocate vectors for GLM statistics.
    p_values = zeros(1, num_features);
    t_values = zeros(1, num_features);
    b_values = zeros(1, num_features);

    % Loop over each FNC feature.
    for f = 1:num_features
        % Extract the f-th FNC feature.
        fnc_feat = fnc(:, f);

        % Combine predictors: the FNC feature and the covariates.
        predictors = [fnc_feat, covariates];

        % Fit a GLM with the specified distribution and link.
        [beta, ~, stats] = glmfit(predictors, clinical_score, dist, 'link', link);

        % Extract statistics for the FNC feature predictor (2nd coefficient, after the intercept).
        p_values(f) = stats.p(2);
        t_values(f) = stats.t(2);
        b_values(f) = beta(2);
    end

    % Apply Benjamini-Hochberg FDR correction.
    q_values = mafdr(p_values, 'BHFDR', true);
end
