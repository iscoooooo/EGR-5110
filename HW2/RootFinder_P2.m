% Written by: Francisco Sanudo
% Date: 2/22/24

function Q = RootFinder_P2(dens, visc, rough, tol)
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


% insert code here

end

%------------------------------------------------------------------------

function [X] = modifiedSecant(FUN,X0,delta,tol)
%{
PURPOSE
modifiedSecant solves systems of nonlinear equations of several variables.

   modifiedSecant attempts to solve equations of the form:

   F(X) = 0    where F and X may be vectors or matrices.

INPUTS
FUN   - set of non-linear algebraic equations of the form F(X) = 0,
        FUN accepts input X and returns a vector (matrix) of equation 
        values F evaluated at X.
X0    - vector of initial guesses
delta - small perturbation
tol   - tolerance

OUPUTS
X - Solution to non-linear equations
%}

% while error > tol

    % Determine value of the function F:
    
    
    
    % Compute the Euclidean norm (relative error):
    
    
    
    % Determine elements of Jacobian J:
    
    
    
    % Solve the linear system -F=JdX for dX at current step
    
    
    
    % Update the value of X


% end

end