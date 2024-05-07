clear; clc; close all;

xi = -10; yi = 5; % initial guesses
f = @(x,y) -10*(x-2).^2 - 5*(y+3).^2 + 20;
tol = 1e-2;
sigma = 0.0001;
beta = 0.5;
% Alter the above lines of code. The function will have only
% one maximum.

% Solver
[xypos,numsteps,numfneval] = calcMaxStudent(f,xi,yi,tol,sigma,beta);

disp('Results:')
disp(' x y')
disp(xypos)
disp(['Number of steps = ', num2str(numsteps)])
disp(['Number of function evaluations = ', num2str(numfneval)])
