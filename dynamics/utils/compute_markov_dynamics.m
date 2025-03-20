function [state_vector, spectral_gap, tv_distance, iterations] = compute_markov_dynamics(m_chain, varargin)
    % compute_markov_dynamics Computes the stationary distribution, spectral gap,
    % and total variation distance convergence of a Markov chain.
    %
    % Usage:
    %   [state_vector, spectral_gap, tv_distance, iterations] = compute_markov_dynamics(m_chain)
    %   [state_vector, spectral_gap, tv_distance, iterations] = compute_markov_dynamics(m_chain, 'mu_initial', mu_init, 'tol', 1e-8, 'max_iter', 1000, 'agg_state_vector', agg_vec)
    %
    % Inputs:
    %   m_chain    - Transition matrix of the Markov chain (n x n)
    %
    % Optional Name-Value Pair Arguments:
    %   'mu_initial'       - Initial distribution vector (n x 1). Defaults to uniform distribution.
    %   'tol'              - Tolerance for convergence. Defaults to 1e-8.
    %   'max_iter'         - Maximum number of iterations. Defaults to 1000.
    %   'agg_state_vector' - Aggregated state distribution vector (n x 1) to use for convergence. Defaults to the computed stationary distribution.
    %
    % Outputs:
    %   state_vector - Stationary distribution vector (n x 1)
    %   spectral_gap - Spectral gap of the transition matrix
    %   tv_distance  - Vector of total variation distances at each iteration
    %   iterations   - Number of iterations taken to converge

%% Parse Optional Inputs
    p = inputParser;
    addParameter(p, 'mu_initial', [], @(x) isnumeric(x) && isvector(x));
    addParameter(p, 'tol', 1e-8, @(x) isnumeric(x) && x > 0);
    addParameter(p, 'max_iter', 10000, @(x) isnumeric(x) && x > 0 && floor(x) == x);
    addParameter(p, 'agg_state_vector', [], @(x) isnumeric(x) && isvector(x));
    parse(p, varargin{:});
    
    mu_initial = p.Results.mu_initial;
    tol = p.Results.tol;
    max_iter = p.Results.max_iter;
    agg_state_vector = p.Results.agg_state_vector;
    
    %% Validate Transition Matrix
    if ~ismatrix(m_chain) || size(m_chain,1) ~= size(m_chain,2)
        error('Transition matrix m_chain must be a square matrix.');
    end
    
    n = size(m_chain, 1);
    
    % Check if rows sum to 1 (within a tolerance)
    row_sums = sum(m_chain, 2);
    if any(abs(row_sums - 1) > 1e-10)
        error('Each row of the transition matrix m_chain must sum to 1.');
    end
    
    % Check for non-negative entries
    if any(m_chain(:) < 0)
        error('Transition matrix m_chain must have non-negative entries.');
    end
    %% Compute Stationary Distribution
    state_vector = compute_stationary_dist(m_chain);
    
    %% If agg_state_vector is provided, validate it
    if ~isempty(agg_state_vector)
        agg_state_vector = agg_state_vector(:); % Ensure it's a column vector
        if length(agg_state_vector) ~= n
            error('agg_state_vector must have the same size as the state space of m_chain.');
        end
        if abs(sum(agg_state_vector) - 1) > 1e-10
            error('agg_state_vector must sum to 1.');
        end
    else
        agg_state_vector = state_vector; % Default to the computed stationary vector
    end

    %% Compute Spectral Gap
    % Compute eigenvalues of P
    eigenvalues = eig(m_chain);
    
    % Sort eigenvalues by their absolute values in descending order
    sorted_eigenvalues = sort(abs(eigenvalues), 'descend');
    
    % The largest eigenvalue should be 1
    lambda_1 = sorted_eigenvalues(1);
    
    if abs(lambda_1 - 1) > 1e-6
        warning('The largest eigenvalue of the transition matrix is not close to 1.');
    end
    
    % Handle cases where there might be multiple eigenvalues equal to 1
    % Find the second distinct eigenvalue
    unique_sorted_eigenvalues = unique(sorted_eigenvalues, 'stable');
    
    if length(unique_sorted_eigenvalues) < 2
        spectral_gap = 0;
        warning('The transition matrix has a spectral gap of 0 (possibly periodic or reducible).');
    else
        lambda_2 = unique_sorted_eigenvalues(2);
        spectral_gap = 1 - lambda_2;
    end
    %% Initialize for Total Variation Distance Computation
    if isempty(mu_initial)
        % Default to uniform distribution
        mu_initial = ones(n,1) / n;
    else
        % Ensure mu_initial is a column vector and sums to 1
        mu_initial = mu_initial(:);
        if abs(sum(mu_initial) - 1) > 1e-10
            mu_initial = mu_initial / sum(mu_initial);
            warning('Initial distribution mu_initial was normalized to sum to 1.');
        end
    end
    
    % Initialize variables
    mu_current = mu_initial;
    tv_distance = zeros(max_iter, 1);
    iterations = 0;
    converged = false;
    
    %% Iterative Process to Compute Total Variation Distance
    for t = 1:max_iter
        % Compute TV distance to the target stationary distribution
        tv_distance(t) = 0.5 * sum(abs(mu_current - agg_state_vector));
        
        % Check for convergence
        if tv_distance(t) < tol
            iterations = t;
            tv_distance = tv_distance(1:t);
            converged = true;
%             fprintf('Converged at iteration %d with TV distance %.2e\n', t, tv_distance(t));
            break;
        end
        
        % Update the current distribution
        mu_current = m_chain' * mu_current;
    end
    
    if ~converged
        iterations = max_iter;
        tv_distance = tv_distance(1:iterations);
        warning('Did not converge within the maximum number of iterations (%d).', max_iter);
    end
end
