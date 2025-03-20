function Hrate = markovEntropyRate(piVec, P)
% markovEntropyRate - Computes the entropy rate of an ergodic Markov chain.
%
% Syntax:
%   Hrate = markovEntropyRate(piVec, P)
%
% Inputs:
%   piVec : 1 x n row vector
%           The stationary distribution of the Markov chain. It should sum to 1 and have non-negative entries.
%
%   P     : n x n matrix
%           The transition probability matrix of the Markov chain, where P(i,j) is the probability
%           of transitioning from state i to state j. Each row must sum to 1, and all entries must be non-negative.
%
% Output:
%   Hrate : Scalar
%           The entropy rate of the Markov chain in bits.
%
% Description:
%   The entropy rate H of an ergodic Markov chain is given by:
%
%       H = - sum_{i=1}^{n} pi(i) * sum_{j=1}^{n} P(i,j) * log2(P(i,j))
%
%   where 0*log2(0) is defined as 0. This function handles the 0*log2(0) case by replacing
%   zero entries in P with ones for the purpose of computing the logarithm.
%
% Example:
%   % Define a 3-state Markov chain.
%   P = [0.7 0.2 0.1;
%        0.3 0.4 0.3;
%        0.2 0.3 0.5];
%
%   % Compute the stationary distribution (using compute_stationary_dist, for example).
%   piVec = compute_stationary_dist(P)';
%
%   % Compute the entropy rate.
%   Hrate = markovEntropyRate(piVec, P);
%
% Author: [Sir-Lord]
% Date: [Today's Date]

    % Check that dimensions of piVec and P are consistent.
    n = length(piVec);
    if size(P,1) ~= n || size(P,2) ~= n
        error('Dimension mismatch: piVec must be length n and P must be n x n.');
    end

    % Validate the stationary distribution piVec.
    if abs(sum(piVec) - 1) > 1e-10
        error('Stationary distribution piVec must sum to 1.');
    end
    if any(piVec < 0)
        error('All elements of the stationary distribution piVec must be non-negative.');
    end

    % Validate the transition matrix P.
    if any(P(:) < 0)
        error('All elements of the transition matrix P must be non-negative.');
    end
    rowSums = sum(P, 2);
    if any(abs(rowSums - 1) > 1e-10)
        error('Each row of the transition matrix P must sum to 1.');
    end

    % Replace zeros in P with ones for the purpose of taking logarithms.
    % This replacement is valid because 0*log2(0) is defined as 0.
    P_nonzero = P;
    P_nonzero(P_nonzero == 0) = 1;

    % Compute the entropy rate.
    Hrate = -sum(piVec .* sum(P .* log2(P_nonzero), 2));
end
