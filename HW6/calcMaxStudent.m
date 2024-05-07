% Written by: Francisco Sanudo
% Date: 5/1/24
%
% PURPOSE
% calcMaxStudent solves a 2D unconstrained optimization problem by finding
% the local maximum of a 2D function f(x,y) using the gradient ascent with
% inexact (backtracking) line search method. The Armijo condition is used
% to ensure the function increases by a minimum amount at each step.
%
% INPUTS
% - f       : Objective function (anonymous func of x and y)
% - xi & yi : Initial guesses for x and y
% - tol     : Error tolerance for convergence
% - sigma   : Armmijo condition constant
% - beta    : Backtracking constant
%
% OUTPUTS
% - xypos    : Array that contains the (x,y) coordinates at each step
% - numsteps : Number of steps required to achieve the termination criteria
% - numfeval : Number of function evaluations
%   --> Everytime f(x,y) is invoked, a counter variable increases by 1
%
% EXAMPLE
% xi = -10; yi = 5;                           % initial guesses
% f = @(x,y) -10*(x-2).^2 - 5*(y+3).^2 + 20;  % objective function
% tol = 1e-2;                                 % tolerance
% sigma = 0.0001;                             % armijo condition consant
% beta = 0.5;                                 % backtracking constant
% 
% [xypos,numsteps,numfneval] = calcMaxStudent(f,xi,yi,tol,sigma,beta);
%
% OTHER
% .m files required              : MAIN.m (calling script)
% Files required (not .m)        : none
% Built-in MATLAB functions used : numel, zeros
% Utility functions              : grad, eucildeanNorm, applyFigureProperties
%
% REFERENCES
% Optimization (notes), P. Nissenson
%

function [xypos,numsteps,numfneval] = calcMaxStudent(f,xi,yi,tol,sigma,beta)

% Initialize Variables
numsteps  = 1;              % Set the iteration counter
numfneval = 0;              % Set the funcion evaluations counter
maxiter   = 10000;          % Maximum iterations
converged = false;          % Convergence condition

X = [xi; yi];               % Set initial condition
xypos = X';                 % Initialize position history

%% Gradient ascent with inexact (backtracking) line search method

% Iterate
while ~converged && numsteps < maxiter
    % Compute gradient
    g = grad(f, X);
    numfneval = numfneval + 4; % Increment function evaluation counter

    % Compute error using gradient norm
    err = euclideanNorm(g);

    if err < tol
        converged = true;
    else
        % Backtracking line search
        h = 1;  % Initial step size

        % Temporary step
        xtemp   = X(1) + g(1)*h;
        ytemp   = X(2) + g(2)*h;

        % Function value at current step
        fcurrent  = f(X(1),X(2));
        numfneval = numfneval + 1; % Increment function evaluation counter

        % Armijo Condition
        while f(xtemp, ytemp) - fcurrent <  sigma*(g'*g)*h
            % Reduce step size using backtracking
            h = beta * h;

            % Increment function evaluation counter
            numfneval = numfneval + 1;
            
            % Update temporary variables
            xtemp = X(1) + g(1)*h;
            ytemp = X(2) + g(2)*h;
        end

        % Update position
        X = X + h*g;

        % Update position history
        xypos = [xypos; X'];

        % Increment iteration counter
        numsteps = numsteps + 1;
    end
end

% Increment function evaluation counter (to account for final step when
% convergence is reached)
numfneval = numfneval + 1;

%% Contour plot that shows the path taken to the maximum

% Define grid
x = -10:.1:10;
y = x;

% Initialize variable
z  = zeros(length(x));

% Evaluate the function at each point on grid
%   Each row corresponds to a constant y value and each column corresponds
%   to a constant x value. This means that z(i, j) represents the function
%   value at x(j) and y(i)
for i = 1:numel(x)
    for j = 1:numel(x)
        z(i,j) = f(x(j), y(i));
    end
end

% Create figure and apply figure properties
f = figure;
position = [0.2, 0.2, 0.5, 0.6];
applyFigureProperties(f, position)

hold on;  % plot over same axes
ax = gca; % get current axes

% Plot contours
contour(x,y,z,'ShowText','on');

% Plot starting point
plot(xi,yi,'bo','Markersize',12,'MarkerFaceColor','b')

% Plot the path to the optimal solution
plot(xypos(:,1),xypos(:,2),'k')
plot(xypos(:,1),xypos(:,2),'ks','MarkerSize',3,'MarkerFaceColor','k')

% Plot the optimal solution
plot(xypos(end,1),xypos(end,2),'m^','Markersize',12,'MarkerFaceColor','m')

% Axis properties
set(ax,'TickLabelInterpreter','latex')
title('Gradient Ascent w/ Armijo Condition')
xlabel('$x$')
ylabel('$y$')
legend('','$x_i$, $y_i$','Path','','$x^*$, $y^*$')
grid on
axis equal

% Colormap & colorbar properties
colormap(ax,'turbo');
cb = colorbar;
cb.TickLabelInterpreter = 'latex';
title(cb,'$\rm Level$','Interpreter','latex')

end

%-------------------------------------------------------------------------

function g = grad(f,X)
% Estimates the gradient of the function f using a centered finite-
% difference approximation.

n = numel(X);          % Dimension of input vector X
g = zeros(n, 1);       % Initialize gradient vector
delta = 0.0001;        % Perturbation parameter

% Estimate gradient
for i = 1:n
    % Create a copy of X to perturb
    X_perturbed = X; 

    % Perturb the i-th element
    X_perturbed(i) = X_perturbed(i) + delta; 

    % Calculate the forward difference
    f_forward = f(X_perturbed(1),X_perturbed(2));

    % Perturb in the negative direction
    X_perturbed(i) = X_perturbed(i) - 2 * delta; 

    % Calculate the backward difference
    f_backward = f(X_perturbed(1),X_perturbed(2));

    % Estimate partial derivative w.r.t. i-th element using central 
    % difference formula
    g(i) = (f_forward - f_backward) / (2 * delta);
end

end

%-------------------------------------------------------------------------

function euclidean_norm_value = euclideanNorm(vec)
% Computes the Euclidean norm of an n-dimensonal vector

% Square each element of vector
squared_elements = vec.^2;

% Initialize sum_of_squares
sum_of_squares = 0;

% Loop through each element
for i = 1:numel(vec)
    % Square the current element
    squared_element = vec(i)^2;

    % Add the squared element to the sum_of_squares
    sum_of_squares = sum_of_squares + squared_element;
end

% Compute the Euclidean norm
euclidean_norm_value = sqrt(sum_of_squares);

end

%-------------------------------------------------------------------------

function applyFigureProperties(figHandle, position)
set(figHandle, ...
    'Units', 'normalized', ...
    'Position', position, ...
    'DefaultTextInterpreter', 'latex', ...
    'DefaultLegendInterpreter', 'latex', ...
    'DefaultAxesFontSize', 14);
end