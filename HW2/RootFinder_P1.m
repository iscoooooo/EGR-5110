% -----------------------MATLAB Function Information-----------------------
%{
Written by: Francisco Sanudo
Date: 2/19/24

PURPOSE
This function solves the non-linear Colebrook equation to determine the
friction factor using the bisection, Newton-Raphson, and secant method.

REFERENCES
Solving Sets of Non-Linear Algebraic Equations (notes), P. Nissenson

INPUTS
- y        : anonymous function *(Colebrook eq.)
- dy       : derivative of the anonymous function y
- fra, frb : lower and upper bound for the Bisection method
- fr0      : initial guess for the Newton-Raphson method
- fr1, fr2 : first and second guess for the Secant method
- tol      : tolerance for termination criteria

OUTPUTS
- frBis   : value of the root using the Bisection method
- iterBis : number of iterations required by the Bisection method
- frNR    : value of the root using the Newton-Raphson method
- iterNR  : number of iterations required by the Newton-Raphson method
- frSec   : value of the root using the Secant method
- iterSec : number of iterations required by the Secant method


OTHER
.m files required              : none
Files required (not .m)        : none
Built-in MATLAB functions used : size, numel, zeros, abs
User-defined functions         : Bisection, newtonRaphson, Secant
%}

function [frBis, iterBis, frNR, iterNR, frSec, iterSec] = ...
    RootFinder_P1(y, dy, fra, frb, fr0, fr1, fr2, tol)

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

set(gcf,'units','normalized','position', [0, 0, .5, .5], ...
    'DefaultTextInterpreter','Latex');
movegui(gcf,'center')

semilogx(tol_values, iter_values, ...
    'LineStyle', ":","Color","blue",'Marker',".",'MarkerSize',12);

title('\textbf{Bisection Method}: Iterations vs Tolerance','FontSize',15);
xlabel('Tolerance','FontSize',15);
ylabel('Iterations','FontSize',15);
grid on;

end

%% Root-Solving Methods

%--------------------------------------------------------------------------

function [X,iter] = Bisection(y,a,b,tol)
i    = 1;   % initialize counter
err  = 1;   % initialize error

% Define search interval
xa   = a; % lower end
xb   = b; % upper end

% Begin loop
while err > tol
    xmid(i) = (xa + xb)/2;

    % determine new bounds
    if y(xa)*y(xmid(i)) > 0
        xa = xmid(i);
    else
        xb = xmid(i);
    end

    xmid(i+1) = (xa + xb)/2;

    err = abs((xmid(i+1) - xmid(i)) / xmid(i+1)); % compute error
    i = i + 1;                                    % increment counter
end

X    = xmid(end); % Final root estimate
iter = i;         % total iterations
end

%--------------------------------------------------------------------------

function [X,iter] = newtonRaphson(y,dy,X0,tol)
i    = 1;   % initialize counter
err  = 1;   % initialize error
x(i) = X0;  % set initial guess

% begin loop
while err > tol
    x(i+1) = x(i) - y(x(i))/dy(x(i));       % compute new root estimate
    err    = abs((x(i+1) - x(i)) / x(i+1)); % compute the eror
    i      = i + 1;                         % increment counter
end

X    = x(end);  % final root estimate
iter = i-1;     % total iterations

end

%--------------------------------------------------------------------------

function [X,iter] = Secant(y,X1,X2,tol)
i    = 2;   % initialize counter
err  = 1;   % initialize error

% set initial guesses
x(i-1) = X1;
x(i)   = X2;

% begin loop
while err > tol
    x(i+1) = x(i) - y(x(i))*(x(i)-x(i-1)) ... % compute new root estimate
             / (y(x(i))-y(x(i-1))); 
    err    = abs((x(i+1) - x(i)) / x(i+1));   % compute the error
    i      = i + 1;                           % increment counter
end

X    = x(end); % final root estimate
iter = i-2;    % total iterations

end