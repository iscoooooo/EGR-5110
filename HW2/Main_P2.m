clc; clear;

% The user inputs the roughness of the pipes, kinematic viscosity of the
% fluid, and termination criterion. I will use different values when 
% testing your code.
dens = 999;         % fluid density, kg/m^3
visc = 1e-3;        % fluid viscosity, N s/m^2
rough = 0.045/1000; % absolute roughness of commercial steel pipe, m
tol = 1E-5;         % tolerance

% The flowrates in RootFinder_P2.m should be a horizontal vector when output.
tic
Q = RootFinder_P2(dens, visc, rough, tol);
toc

% The transpose operator converts Qstudent into a column vector.
Q = Q';