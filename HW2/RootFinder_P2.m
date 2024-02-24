
% Written by: Francisco Sanudo
% Date: 2/22/24

function Q = RootFinder_P2()
%{
PURPOSE
This function calculates the volumetric flowrates through each branch of a 
pipe network using the Modified Secant method -- a method that uses a
finite-differencing approximation of first-order derivatives to solve a 
set of n non-linear algebraic equations with n unknowns.

REFERENCES
Solving Sets of Non-Linear Algebraic Equations (notes), P. Nissenson

INPUTS
- Q : An array contaning the flowrates in each branch

OUTPUTS
- dens  : density of the fluid
- visc  : kinematic viscosity of the fluid
- rough : absolute pipe roughness
- tol   : tolerance for termination criteria


OTHER
.m files required              : none
Files required (not .m)        : none
Built-in MATLAB functions used : size, numel, zeros, abs
User-defined functions         : 
%}

X0    = [1;1];
F     = @(X) [X(1)^2 + X(1)*X(2) - 10; X(2) + 3*X(1)*X(2)^2 - 57];
delta = 0.01;
tol   = 1e-5;

Q = multivariateRoot(F,X0,delta,tol);

end

%------------------------------------------------------------------------

function [X] = multivariateRoot(FUN,X0,delta,tol)
%{
PURPOSE
multivariate solves systems of nonlinear equations of several variables.

   multivariate attempts to solve equations of the form:

   F(X) = 0    where F and X may be vectors or matrices.

INPUTS
FUN   - set of non-linear algebraic equations of the form F(X) = 0,
        FUN accepts input X and returns a vector (matrix) of equation 
        values F evaluated at X.
X0    - vector of initial guesses
delta - small perturbation
tol   - tolerance

OUPUTS
X - Solution to set of non-linear equations
%}

% Pre-processing
k     = 1;
error = 1;
X(:,k)  = X0;

% Begin loop
while error > tol

    % Determine current value of the function F:
    F(:,k) = FUN(X(:,k));

    % Compute the Euclidean norm (relative error):
    error = sqrt(F(:,k)'*F(:,k));

    % Determine elements of Jacobian J:
    J = Jacobian(FUN,X(:,k),delta);

    % Solve the linear system -F=JdX for dX at current step
    dX = -J\F(:,k);

    % Update the value of X
    X(:,k+1) = X(:,k) + dX;

    % Increment counter
    k = k + 1;
end

X = X(:,end);

end

function J = Jacobian(FUN,X0,delta)
% This function estimates the Jacobian of a multivariate, vector-valued
% function using Finite Differencing.

F = FUN;
m = numel(X0);      % number of rows
n = numel(F(X0)); % number of columns

% Intialize Jacobian matrix
J = zeros(m,n);

% Begin loop
for k = 1:n

    % Only perturb the k-th component of X0
    X_perturbed = X0;
    X_perturbed(k) = X_perturbed(k) + delta;

    % Estimate the k-th column of the Jacobian
    J(:,k) = (F(X_perturbed) - F(X0)) / delta;

end

end