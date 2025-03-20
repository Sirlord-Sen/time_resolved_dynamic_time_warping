function isErgodic = checkIsErgodic(P)
% checkIsErgodic - Checks if a Markov chain is ergodic.
%
% Syntax:
%   isErgodic = checkIsErgodic(P)
%
% Inputs:
%   P - n x n transition matrix of the Markov chain.
%
% Outputs:
%   isErgodic - Logical value (true/false) indicating if the Markov chain is
%               ergodic (i.e., irreducible and aperiodic).
%
% Description:
%   A Markov chain is ergodic if it is both irreducible and aperiodic. This function
%   first checks that P is a valid transition matrix (square, non-negative, and rows 
%   summing to 1). It then creates a directed graph from P to assess irreducibility
%   by verifying that the graph is strongly connected. If irreducible, the function
%   checks for aperiodicity by examining the return times to a reference state and 
%   computing the greatest common divisor (GCD) of these return times. If the GCD is 1,
%   the chain is aperiodic.
%
% Example:
%   P = [0.7 0.2 0.1;
%        0.3 0.4 0.3;
%        0.2 0.3 0.5];
%   isErgodic = checkIsErgodic(P);
%
% Author: Sir-Lord

    % Validate input: Ensure P is square.
    [n, m] = size(P);
    if n ~= m
        error('Transition matrix P must be square.');
    end

    % Check that each row of P sums to 1 (within tolerance).
    if any(abs(sum(P, 2) - 1) > 1e-10)
        isErgodic = false;
        return;
    end

    % Create a directed graph from the transition matrix (nonzero entries indicate edges).
    G = digraph(P > 0);
    
    % Check for irreducibility: The Markov chain is irreducible if the graph is strongly connected.
    scc = conncomp(G, 'Type', 'strong');
    if max(scc) > 1
        isErgodic = false;
        return;
    end
    
    % Check for aperiodicity.
    % For a finite Markov chain, examine the return times to a reference state.
    max_k = 2 * n;  % Upper bound on steps to check.
    Pi = P;
    state = 1;  % Choose state 1 as the reference state.
    return_steps = [];
    
    for k = 1:max_k
        if Pi(state, state) > 1e-12  % A return is detected.
            return_steps(end+1) = k; %#ok<AGROW>
        end
        Pi = Pi * P;
    end
    
    if isempty(return_steps)
        % No return detected within max_k steps.
        isErgodic = false;
        return;
    end
    
    % Compute the greatest common divisor (GCD) of the return steps.
    d = return_steps(1);
    for idx = 2:length(return_steps)
        d = gcd(d, return_steps(idx));
        if d == 1
            break;  % If GCD is 1, the chain is aperiodic.
        end
    end
    
    % The chain is aperiodic if the GCD of the return times is 1.
    isAperiodic = (d == 1);
    
    % The chain is ergodic if it is both irreducible (already verified) and aperiodic.
    isErgodic = isAperiodic;
end
