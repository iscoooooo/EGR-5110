% Written by: Francisco Sanudo
% Date: 2/22/24
%
% PURPOSE
% RootFinder_P2 calculates the volumetric flowrates through each branch of a 
% pipe network using the Modified Secant method -- a method that uses a
% finite-differencing approximation of first-order derivatives to solve a 
% set of n non-linear algebraic equations with n unknowns.
% 
% REFERENCES
% Solving Sets of Non-Linear Algebraic Equations (notes), P. Nissenson
% 
% INPUTS
% - Q : An array contaning the flowrates in each branch
% 
% OUTPUTS
% - dens  : density of the fluid
% - visc  : kinematic viscosity of the fluid
% - rough : absolute pipe roughness
% - tol   : tolerance for termination criteria
% 
% OTHER
% .m files required              : none
% Files required (not .m)        : none
% Built-in MATLAB functions used : size, numel, zeros, abs
% User-defined functions         : myFunc, Jacobian, multivariateRootFinder

function Q = RootFinder_P2(dens, visc, rough, tol)

% check if inputs are numeric
if ~isnumeric(dens) || ~isnumeric(visc) || ~isnumeric(rough) || ~isnumeric(tol)
    error('All inputs must be numeric.');
end

% check if inputs are scalar
if ~isscalar(dens) || ~isscalar(visc) || ~isscalar(rough) || ~isscalar(tol)
    error('All inputs must be scalar.');
end

% Ensure positive values input parameters
if dens <= 0 || visc <= 0 || rough <= 0 || tol <= 0
    error('Density, viscosity, roughness, and tolerance must be positive.');
end

% Pipe lengths
L = [200;100;360;200;100;200;300;450]; % [meters]

% Pipe diameters
D = [123.4;158.6;123.4;123.4;176.2;96.8;123.4;109.8]/1000; % [meters]

% Initial guess
Q0 = [1200;800;300;900;100;300;300;400]/3600; % [meters^3/sec]

% Perturbation (arbitrary small amount) for estimating Jacobian
delta = 0.01;

% System parameters
parameters = struct('Lengths',L,'Diameters',D,'Density',dens,'Viscosity',visc, ...
    'Roughness',rough);

% Create anonymous function with additional parameters
myEquations = @(Q) myFunc(Q,parameters);

% Determine flowrates
[Q,~,QPath,normPath] = multivariateRootFinder(myEquations,Q0,delta,tol);

% Q(1) at each time step
Q1_path = QPath(1,:);

% Plot Q(1) and norm at each iteration
plotResults(Q1_path, normPath)

end

%% Function handle for set of non-linear equations

function F = myFunc(Q,parameters)
% myFunc Returns the value of the vector-valued function that defines the
% set of non-linear algebraic equations for a pipe network
%
%
% Inputs:
% - Q: An array contaning pipe flowrates at each iteration
% - parameters: A structure with fields containing various system parameters

% Parameters
L     = parameters.Lengths;   % [m]
D     = parameters.Diameters; % [m]
dens  = parameters.Density;   % [kg/m^3]
visc  = parameters.Viscosity;      % [N*s/m^2]
rough = parameters.Roughness; % [m]
g     = 9.81;                 % [m/s^2]

% Calculate Reynold's number for each pipe
Re = 4.*abs(Q)*dens./(pi.*D.*visc);

% Calculate friction factor for each pipe (Haaland Eq.)
fr = (-1.8.*log10(((rough./D)./3.7).^1.11 + 6.9./Re)).^-2;

% Set of non-linear algebraic equations
F = [Q(1) + Q(2) - 2000/3600; ...                       % F(1)

Q(3) + Q(7) + Q(8) - 1000/3600; ...                % F(2)

-Q(4) - Q(5) - Q(7) + 1300/3600; ...                % F(3)

-Q(1) + Q(5) + Q(6) + 800/3600; ...                 % F(4)

-Q(6) - Q(8) + 700/3600; ...                        % F(5)

-(fr(1)*8*L(1)/(pi^2*g*D(1)^5))*Q(1)*abs(Q(1)) + ... % F(6)
(fr(2)*8*L(2)/(pi^2*g*D(2)^5))*Q(2)*abs(Q(2)) + ...
(fr(4)*8*L(4)/(pi^2*g*D(4)^5))*Q(4)*abs(Q(4)) - ...
(fr(5)*8*L(5)/(pi^2*g*D(5)^5))*Q(5)*abs(Q(5));  ...

-(fr(3)*8*L(3)/(pi^2*g*D(3)^5))*Q(3)*abs(Q(3)) - ... % F(7)
(fr(4)*8*L(4)/(pi^2*g*D(4)^5))*Q(4)*abs(Q(4)) + ...
(fr(7)*8*L(7)/(pi^2*g*D(7)^5))*Q(7)*abs(Q(7));  ...

(fr(5)*8*L(5)/(pi^2*g*D(5)^5))*Q(5)*abs(Q(5)) - ... % F(8)
(fr(6)*8*L(6)/(pi^2*g*D(6)^5))*Q(6)*abs(Q(6)) - ...
(fr(7)*8*L(7)/(pi^2*g*D(7)^5))*Q(7)*abs(Q(7)) + ...
(fr(8)*8*L(8)/(pi^2*g*D(8)^5))*Q(8)*abs(Q(8))];

end

%% Main Algorithm

function [X,iter,XPath,normPath] = multivariateRootFinder(FUN,X0,delta,tol)
% PURPOSE
% multivariateRootFinder solves systems of nonlinear equations of several variables.
% 
%    multivariateRootFinder attempts to solve equations of the form:
% 
%    F(X) = 0    where F and X may be vectors or matrices.
% 
% INPUTS
% FUN   - set of non-linear algebraic equations of the form F(X) = 0,
%         FUN accepts input X and returns a vector (matrix) of equation 
%         values F evaluated at X.
% X0    - vector of initial guesses
% delta - small perturbation
% tol   - tolerance
% 
% OUPUTS
% X - Solution to set of non-linear equations

% Input validation
if ~isa(FUN, 'function_handle')
    error('multivariateRootFinder:InvalidInput', 'FUN must be a function handle.');
end

if ~isnumeric(X0) || ~isvector(X0)
    error('multivariateRootFinder:InvalidInput', 'X0 must be a numeric vector.');
end

if ~isnumeric(delta) || numel(delta) ~= 1 || delta == 0
    error('multivariateRootFinder:InvalidInput', 'delta must be a non-zero numeric scalar.');
end

try
    F0 = FUN(X0);
catch ME
    error('multivariateRootFinder:FunctionEvaluationError', 'Error evaluating FUN at X0: %s', ME.message);
end

if ~isnumeric(F0) || ~isvector(F0)
    error('multivariateRootFinder:InvalidFunctionOutput', 'FUN must return a numeric vector.');
end

% Initialize variables
max_iterations = 1000;  % maximum iterations
converged      = false; % convergence condition
norm           = 1;     % initialize norm
k              = 0;     % set counter
X              = X0;    % initial condition

% Begin loop
while k < max_iterations && ~converged

    % Increment counter
    k = k + 1;

    % Determine current value of the function F:
    F = FUN(X(:,k));

    % Compute the Euclidean norm (relative error):
    norm(k+1) = sqrt(F'*F);

    if norm(k+1) < tol
        converged = true;
    end

    % Determine elements of Jacobian J:
    J = Jacobian(FUN,X(:,k),delta);

    % Solve the linear system -F=JdX for dX at current step
    dX = -J\F;

    % Update the value of X
    X(:,k+1) = X(:,k) + dX;

    % Check for numerical instability
    if any(isnan(X(:,k+1))) || any(isinf(X(:,k+1)))
        error('multivariateRootFinder:NumericalInstability', 'The computation resulted in NaN or Inf.');
    end

end

% check for convergence
if ~converged
    warning('multivariateRootFinder:NonConvergence', 'The function did not converge to a solution.');
end

XPath    = X(:,:);    % Entire history of X at each step
normPath = norm(:);   % Entire history of norm at each step
iter     = k-1;       % Total iterations
X        = X(:,end)'; % Solution

end

%% Jacobian estimate function

function J = Jacobian(FUN,X0,delta)
% This function estimates the Jacobian of a multivariate, vector-valued
% function using Finite Differencing.
%
% Input validation is ignored since the function multivariateRootFinder
% performs the necessary check.

% Initialize variables
F = FUN;
m = numel(X0);    % number of rows
n = numel(F(X0)); % number of columns

% Intialize Jacobian matrix
J = zeros(m,n);

% Begin loop
for k = 1:n

    % Only perturb the k-th component of X0
    X_perturbed = X0;
    X_perturbed(k) = X_perturbed(k) + delta;

    % Error handling
    try
        F_perturbed = FUN(X_perturbed);
    catch ME
        error('Jacobian:FunctionEvaluationError', 'Error evaluating FUN at X_perturbed: %s', ME.message);
    end

    % Ensure that FUN returns consistent sizes
    if numel(F_perturbed) ~= n
        error('Jacobian:InconsistentOutput', 'FUN must return outputs of consistent size.');
    end

    % Estimate the k-th column of the Jacobian
    J(:,k) = (F(X_perturbed) - F(X0)) / delta;

end

end

%% Plotting function

function plotResults(Q1_values, norms)
% plotIterationResults Plots the Q(1) values and Euclidean norms across iterations
%
% Inputs:
% - Q1_values: An array of Q(1) values at each iteration
% - norms: An array of Euclidean norms at each iteration

figure;
% figure settings
set(gcf,'units','normalized','position', [0, 0, .4, .5], ...
    'DefaultTextInterpreter','Latex');
movegui(gcf,'center')

% First subplot for Q(1)
subplot(2, 1, 1);
hold on
plot(Q1_values);
yline(Q1_values(end),'k--')
hold off
title('\textbf{$Q_1$ vs. Iteration Number}','FontSize',18);
xlabel('Iteration','FontSize',18);
ylabel('$Q_1 \,\mathrm{\frac{m^3}{s}}$','FontSize',18);
ylim([0, max(Q1_values) + .10*max(Q1_values)])
legend('$Q_1$ path','Solution','interpreter','latex')

% Second subplot for the Euclidean norm
subplot(2, 1, 2);
plot(norms);
title('\textbf{Euclidean Norm vs. Iteration Number}','FontSize',18);
xlabel('Iteration','FontSize',18);
ylabel('$\|\mathbf{F}\|_e$','FontSize',18);

end
