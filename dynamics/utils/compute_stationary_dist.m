function state_vector = compute_stationary_dist(m_chain)
% compute_stationary_dist - Computes the stationary distribution of a Markov chain.
%
% Syntax:
%   state_vector = compute_stationary_dist(m_chain)
%
% Inputs:
%   m_chain - A square transition probability matrix (n x n) of the Markov chain.
%
% Outputs:
%   state_vector - A column vector (n x 1) representing the stationary distribution,
%                  i.e., the unique probability vector pi satisfying:
%                        pi' * m_chain = pi'  and  sum(pi) = 1.
%
% Description:
%   This function computes the stationary distribution of a Markov chain by solving
%   the linear system:
%
%       pi' * m_chain = pi'
%
%   subject to the normalization condition:
%
%       sum(pi) = 1.
%
%   The system is reformulated as:
%
%       (m_chain' - I) * pi = 0,
%
%   and then the normalization condition is appended. The resulting system is solved
%   using MATLAB's backslash operator.
%
% Example:
%   % Define a transition matrix
%   P = [0.7 0.2 0.1;
%        0.3 0.4 0.3;
%        0.2 0.3 0.5];
%
%   % Compute the stationary distribution
%   pi = compute_stationary_dist(P);
%
% Author: [Sir-Lord]

    n = size(m_chain, 1);
    % Transpose the matrix to set up the equation pi' * P = pi'
    PT = m_chain';
    % Form the homogeneous system: (P' - I) * pi = 0
    A = PT - eye(n);
    
    % Append the normalization condition sum(pi) = 1
    A = [A; ones(1, n)];
    b = [zeros(n, 1); 1];
    
    % Solve the linear system
    pi_solution = A \ b;
    
    % Ensure the solution is real (remove any negligible imaginary parts)
    pi_solution = real(pi_solution);
    
    % Normalize to ensure the solution sums to 1
    state_vector = pi_solution / sum(pi_solution);
end
