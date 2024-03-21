function [X,iter] = Secant(y,X1,X2,tol)
max_iterations = 1000;  % maximum iterations
converged      = false; % convergence condition
i = 1;                  % initialize counter

% set initial guesses
x(i) = X1;
x(i+1) = X2;

% begin loop
while i < max_iterations && ~converged
    i = i + 1; % increment counter

    % compute new root estimate
    x(i+1) = x(i) - y(x(i))*(x(i)-x(i-1)) ... 
        / (y(x(i))-y(x(i-1)));

    % compute the error
    err = abs((x(i+1) - x(i)) / x(i+1));
    if err < tol
        converged = true;
    end

    % check for numerical instability
    if any(isnan(x(i+1))) || any(isinf(x(i+1)))
        error('Secant:NumericalInstability', ['The computation resulted ' ...
            'in NaN or Inf.']);
    end

end

% check for convergence
if ~converged
    warning('Secant:NonConvergence', ['The function did not ' ...
        'converge to a solution.']);
end

X    = x(end); % final root estimate
iter = i-1;    % total iterations

end