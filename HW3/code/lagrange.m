clc;clear;close all

% Define constants and linspace
mu = 1/82.45;
x = linspace(-4,4,10000);

% Calculate the function values to visualize possible roots
y1 = x - (1-mu)./(x+mu).^2 + mu./(x-1+mu).^2; % L1
y2 = x - (1-mu)./(x+mu).^2 - mu./(x-1+mu).^2; % L2
invalid = x + (1-mu)./(x+mu).^2 - mu./(x-1+mu).^2; % invalid region
y3 = x + (1-mu)./(x+mu).^2 + mu./(x-1+mu).^2; % L3

% Create figure
fig = figure;
applyFigureProperties(fig,[0.2,0.2,0.5,0.6])

hold on
yline(0,'w--')
h1 = plot(x,y1,'r', 'LineWidth',2, 'DisplayName','L1');
h2 = plot(x,y2,'g', 'LineWidth',2, 'DisplayName','L2');
h3 = plot(x,invalid,'y', 'LineWidth',2, 'DisplayName','Invalid');
h4 = plot(x,y3,'b', 'LineWidth',2, 'DisplayName','L3');
h5 = plot(-mu,0,'co', ...
    'MarkerSize',12, ...
    'MarkerFaceColor','c', ...
    'MarkerEdgeColor','c',...
    'DisplayName','Earth');
h6 = plot(1-mu,0,'wo', ...
    'MarkerSize',12, ...
    'MarkerFaceColor','w', ...
    'MarkerEdgeColor','w',...
    'DisplayName','Moon');

% Set graph properties
title('Co-linear Solutions')
legend([h1,h2,h3,h4,h5,h6],'Location','northwest')
grid on
xlim([-2,2]), ylim([-4, 4])

applyAxisAndLegendProperties(fig)

%% Solve for Collinear Lagrange Points

% redefine this as a set of non-linear equations and use fsolve

% Define anyonymous functions
func1 = @(x) x - (1-mu)/(x+mu)^2 + mu/(x-1+mu)^2; % L1 equation
func2 = @(x) x - (1-mu)/(x+mu)^2 - mu/(x-1+mu)^2; % L2 equation
func3 = @(x) x + (1-mu)/(x+mu)^2 + mu/(x-1+mu)^2; % L3 equation

% Define intial guesses
x1_L1 = 0.6;  x2_L1 = 0.7;
x1_L2 = 1.0;  x2_L2 = 1.1;
x1_L3 = -0.5; x2_L3 = -0.4;

% Set tolerance
tol = 1e-12;

% Use secant method to find Lagrange points
L1 = Secant(func1,x1_L1,x2_L1,tol);
L2 = Secant(func2,x1_L2,x2_L2,tol);
L3 = Secant(func3,x1_L3,x2_L3,tol);

% Plot lagrange points on existing figure
plot(L1,0,'ro','MarkerFaceColor','r','MarkerSize',10,'HandleVisibility','off')
plot(L2,0,'go','MarkerFaceColor','g','MarkerSize',10,'HandleVisibility','off')
plot(L3,0,'bo','MarkerFaceColor','b','MarkerSize',10,'HandleVisibility','off')

%% Solve for Non-Collinear Lagrange Points

xguess = 0.8;

A = [-(1-mu)*(xguess + mu), -mu*(xguess-1+mu);
    mu-1, -mu];
b = [xguess;1];

sol = -A^-1*b;

alpha = sol(1);
beta  = sol(2);