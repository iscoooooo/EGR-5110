% Written by: Francisco Sanudo
% Date: 2/19/24
%
% PURPOSE
% This function solves the non-linear Colebrook equation to determine the
% friction factor using the bisection, Newton-Raphson, and secant method.
%
% REFERENCES
% Solving Sets of Non-Linear Algebraic Equations (notes), P. Nissenson
%
% INPUTS
% - y        : anonymous function *(Colebrook eq.)
% - dy       : derivative of the anonymous function y
% - fra, frb : lower and upper bound for the Bisection method
% - fr0      : initial guess for the Newton-Raphson method
% - fr1, fr2 : first and second guess for the Secant method
% - tol      : tolerance for termination criteria
%
% OUTPUTS
% - frBis   : value of the root using the Bisection method
% - iterBis : number of iterations required by the Bisection method
% - frNR    : value of the root using the Newton-Raphson method
% - iterNR  : number of iterations required by the Newton-Raphson method
% - frSec   : value of the root using the Secant method
% - iterSec : number of iterations required by the Secant method
%
%
% OTHER
% .m files required              : none
% Files required (not .m)        : none
% Built-in MATLAB functions used : size, numel, zeros, abs
% User-defined functions         : Bisection, newtonRaphson, Secant

function [frBis, iterBis, frNR, iterNR, frSec, iterSec] = ...
    RootFinder_P1(y, dy, fra, frb, fr0, fr1, fr2, tol)
%% Input Validation

% Check if inputs are numeric
if ~isnumeric(fra) || ~isnumeric(frb) || ~isnumeric(fr0) || ...
        ~isnumeric(fr1) || ~isnumeric(fr2) || ~isnumeric(tol)
    error('RootFinder_P1:InvalidInput', 'Arguments must be a numeric value.');
end

% Check if inputs are scalar
if ~isscalar(fra) || ~isscalar(frb) || ~isscalar(fr0) || ...
        ~isscalar(fr1) || ~isscalar(fr2) || ~isscalar(tol)
    error('RootFinder_P1:InvalidInput', 'Arguments must be a scalar.');
end

% Check if y is a function handle
if isa(y, 'function_handle')
    % Additional check to confirm it's an anonymous function
    info = functions(y);
    if ~isfield(info, 'function') && isempty(info.function)
        error('RootFinder_P1:InvalidInput', ['The input y is a function handle,' ...
            ' but not an anonymous function.']);
    end
else
    error('RootFinder_P1:InvalidInput', 'The input y must be an anonymous function.')
end

% Check if dy is a function handle
if isa(dy, 'function_handle')
    % Additional check to confirm it's an anonymous function
    info = functions(dy);
    if ~isfield(info, 'function') && isempty(info.function)
        error('RootFinder_P1:InvalidInput', ['The input dy is a function handle, ' ...
            'but not an anonymous function.']);
    end
else
    error('RootFinder_P1:InvalidInput', 'The input dy must be an anonymous function.')
end

%% Bisection Method
[frBis,iterBis] = Bisection(y,fra,frb,tol);

%% Newton-Raphson Method
[frNR,iterNR] = newtonRaphson(y,dy,fr0,tol);

%% Secant Method
[frSec,iterSec] = Secant(y,fr1,fr2,tol);

%% Bisection Method Plot
tol_values  = logspace(-10, -5, 100);   % Generate 100 log-spaced tol values 1e-10 and 1e-5
iter_values = zeros(size(tol_values));  % Initialize array to store iter values

for k = 1:numel(tol_values)
    tol_k = tol_values(k);

    % Perform bisection method with the current tol value
    [~,i] = Bisection(y,fra,frb,tol_k);

    iter_values(k) = i;
end

% Plot the results
figure;

set(gcf,'units','normalized','position', [0, 0, .4, .5], ...
    'DefaultTextInterpreter','Latex');
movegui(gcf,'center')

semilogx(tol_values, iter_values, ...
    'LineStyle', ":","Color","blue",'Marker',".",'MarkerSize',12);

title('\textbf{Bisection Method}: Iterations vs Tolerance','FontSize',15);
xlabel('Tolerance','FontSize',15);
ylabel('Iterations','FontSize',15);
legend('Iterations','fontsize',14,'interpreter','latex')
grid on;

end

%% Root-Solving Methods

%--------------------------------------------------------------------------

function [X,iter] = Bisection(y,a,b,tol)
max_iterations = 1000;  % maximum iterations
converged      = false; % convergence condition
i = 0;                  % initialize counter

% Define search interval
xa   = a; % lower end
xb   = b; % upper end

% Begin loop
while i < max_iterations && ~converged
    i = i + 1;             % increment counter

    % current root estimate
    xmid(i) = (xa + xb)/2;

    % determine new bounds
    if y(xa)*y(xmid(i)) > 0
        xa = xmid(i);
    else
        xb = xmid(i);
    end

    % compute next root estimate
    xmid(i+1) = (xa + xb)/2;

    % compute error
    err = abs((xmid(i+1) - xmid(i)) / xmid(i+1));
    if err < tol
        converged = true;
    end

    % check for numerical instability
    if any(isnan(xmid(i+1))) || any(isinf(xmid(i+1)))
        error('Bisection:NumericalInstability', ['The computation resulted ' ...
            'in NaN or Inf.']);
    end
end

% check for convergence
if ~converged
    warning('Bisection:NonConvergence', ['The function did not ' ...
        'converge to a solution.']);
end

X    = xmid(end); % final root estimate
iter = i;         % total iterations

end

%--------------------------------------------------------------------------

function [X,iter] = newtonRaphson(y,dy,X0,tol)
max_iterations = 1000;  % maximum iterations
converged      = false; % convergence condition
x = X0;                 % set initial guess
i = 0;                  % initialize counter

% begin loop
while i < max_iterations && ~converged
    i = i + 1; % increment counter

    % compute new root estimate
    x(i+1) = x(i) - y(x(i))/dy(x(i));
    
    % compute the eror
    err    = abs((x(i+1) - x(i)) / x(i+1)); 
    if err < tol
        converged = true;
    end

    % check for numerical instability
    if any(isnan(x(i+1))) || any(isinf(x(i+1)))
        error('newtonRaphson:NumericalInstability', ['The computation resulted ' ...
            'in NaN or Inf.']);
    end
end

% check for convergence
if ~converged
    warning('newtonRaphson:NonConvergence', ['The function did not ' ...
        'converge to a solution.']);
end

X    = x(end); % final root estimate
iter = i;      % total iterations

end

%--------------------------------------------------------------------------

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