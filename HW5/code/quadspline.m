% Written by: Francisco Sanudo
% Date: 4/26/24
%
% PURPOSE
% quadspline interpolates velocity data using quadratic splines and provide
% computations related to these splines, including:
%
%   1. Interpolation: Fit quadratic splines to given time and velocity data
%   points
%
%   2. Integration: Compute the integral of velocity over specified time
%   intervals to determine distances traveled
%
%   3. Velocity Query: Determine the instantaneous velocity at a queried
%   time
%  
%   4. Distance Calculation: Compute the distance traveled between two
%   specified times
%
% INPUTS
% - x        : Time vector corresponding to velocity data
% - fx       : Velocity vector
% - tinstant : Time at which to determine the instantaneous velocity
% - t1 & t2  : Initial and final times for distance calculations
%
% OUTPUTS
% - totaldist : Total distance traveled from the beginning to the end of the time vector
% - vinstant  : Instantaneous velocity at the specified time 'tinstant'
% - subdist   : Distance traveled between the specified times 't1' and 't2'
%
% EXAMPLE
% x = [0, 1, 2, 3];               % Time vector
% fx = [0, 2, 1, 3];              % Velocity vector
% tinstant = 1.5;                 % Instantaneous velocity query time
% t1 = 0.5;                       % Start time for distance calculation
% t2 = 2.5;                       % End time for distance calculation
%
% [totaldist, vinstant, subdist] = quadspline(x, fx, tinstant, t1, t2);
%
% OTHER
% .m files required              : MAIN.m (calling script)
% Files required (not .m)        : none
% Built-in MATLAB functions used : numel, zeros, fplot
% User-defined functions         : applyFigureProperties
%
% REFERENCES
% Numerical Integration (notes), P. Nissenson
%

function [totaldist,vinstant,subdist] = quadspline(x,fx,tinstant,t1,t2)

% Initialize variables
N = numel(x) - 1;   % total number of splines
b = zeros(3*N,1);   % column matrix of knowns
A = zeros(3*N,3*N); % coefficient matrix

%% Generate coefficients of splines

j = 1; % spline index
k = 1; % starting column index

% Condition #1: Functions are continuous at interior knots
for i = 2:2:2*N
    A(i,k:k+2)   = [x(j)^2, x(j), 1];
    b(i)         = fx(j);
    j            = j + 1;             % increment spline number
    A(i+1,k:k+2) = [x(j)^2, x(j), 1];
    b(i+1)       = fx(j);
    k            = k + 3;             % increment starting column index
end

j = 1; % spline index [reset]
k = 1; % starting column index [reset]

% Condition #2: First and last functions pass through end knots
%   <Satisfied in loop above>

% Condition #3: First derivatives are continuous at interior knots
for i = 2*N+2:3*N
    A(i,k:k+4) = [2*x(j), 1, 0, -2*x(j), -1];
    j          = j + 1;
    k          = k + 3;
end

% Condition #4: Assume second derivative is 0 at first knot
A(1,1) = 1;

% Solve for unknown coefficients
c = A\b;

%% Plot splines

j = 1;  % spline index [reset]

% Create figure and apply figure properties
f = figure;
position = [0.2, 0.2, 0.5, 0.6];
applyFigureProperties(f, position)

hold on % plot over same axes

% Begin loop to plot all splines
for i = 1:N
    spline = @(X) c(j)*X.^2 + c(j+1).*X + c(j+2); % create spline function handle
    fplot(spline,[x(i) x(i+1)])                   % plot current spline between x(i) and x(i+1)

    j = j + 3; % increment to next set of coefficients
end

% Plot data points
h = plot(x,fx,'ro');

% Axis properties
set(gca,'TickLabelInterpreter','latex')
title('Quadratic Spline Fit')
xlabel('Time ($s$)');
ylabel('Distance ($m$)')
legend(h, 'Data Points')
grid on

%% Calculate total distance traveled

% Initialize variables
totaldist = 0;
k = 1;          % coefficient index

% Integration of the velocity function over the current segment

for j = 2:N+1 % N+1 data points where N is the number of splines
    I = (c(k)*x(j)^3/3 + c(k+1)*x(j)^2/2 + c(k+2)*x(j)) ...
        - (c(k)*x(j-1)^3/3 + c(k+1)*x(j-1)^2/2 + c(k+2)*x(j-1));
    totaldist = totaldist + I;

    k = k + 3; % increment to next set of coefficients
end

%% Calculate velocity at tinstant

% Find the relevant spline that contains tinstant
spline_idx = findSpline(x,tinstant);

% Starting coefficient index
k = 3*spline_idx - 2;

% Calculate the velocity at tinstant
vinstant = c(k)*tinstant^2 + c(k+1)*tinstant + c(k+2);

%% Calculate the distance traveled from t1 to t2

% Find the relevant spline that contains t1 and t2
spline_idx_t1 = findSpline(x,t1);
spline_idx_t2 = findSpline(x,t2);

% Initialize variables for distance calculation
subdist = 0;
k = 3*spline_idx_t1 - 2;  % Starting coefficient index for t1 segment

% Loop through spline segments from t1 to t2
for j = (spline_idx_t1):(spline_idx_t2)

    % Integration of the velocity function over the current segment

    % [t1, x(j+1)]
    if j == spline_idx_t1
        I = (c(k)*x(j+1)^3/3 + c(k+1)*x(j+1)^2/2 + c(k+2)*x(j+1)) ...
            - (c(k)*t1^3/3 + c(k+1)*t1^2/2 + c(k+2)*t1);
        subdist = subdist + I;
    % [x(j), t2]
    elseif j == spline_idx_t2
        I = (c(k)*t2^3/3 + c(k+1)*t2^2/2 + c(k+2)*t2) ...
            - (c(k)*x(j)^3/3 + c(k+1)*x(j)^2/2 + c(k+2)*x(j));
        subdist = subdist + I;
    % [x(j), x(j+1)]    
    else 
        I = (c(k)*x(j+1)^3/3 + c(k+1)*x(j+1)^2/2 + c(k+2)*x(j+1)) ...
            - (c(k)*x(j)^3/3 + c(k+1)*x(j)^2/2 + c(k+2)*x(j));
        subdist = subdist + I;
    end

    % Update coefficient index for the next segment
    k = k + 3;
end

end

%------------------------------------------------------------------------

function idx = findSpline(x, tinstant)
% Finds the spline segment index that contains tinstant

% Initialize segment index
idx = 0;

% Check each interval [x(i), x(i+1)] to find the relevant segment
for i = 1:numel(x)-1
    if tinstant >= x(i) && tinstant < x(i+1)
        idx = i;
        break;  % Exit loop once segment is found
    end
end

% If tinstant is out of range (before first element or after last element)
if isempty(idx)
    error('tinstant is out of the range of input time vector x.');
end

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