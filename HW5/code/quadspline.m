% Written by: Francisco Sanudo
% Date: 4/26/24
%
% PURPOSE
% quadspline generates a set of quadratic splines connecting the data
% points and calculates the integral under the splines.
%
% REFERENCES
% Numerical Integration (notes), P. Nissenson
%
% INPUTS
% - x        : time vector
% - fx       : velocity vector
% - tinstant : time queried to determine the instantaneous velocity
% - t1       : initial time
% - t2       : final time
%
% OUTPUTS
% - totaldist : total distance traveled from initial time to final time
% - vinstant  : velocity at tinstant
% - subdist   : distance traveled from t1 to t2
%
% OTHER
% .m files required              : MAIN.m (calling script)
% Files required (not .m)        : none
% Built-in MATLAB functions used : numel, zeros, fplot
% User-defined functions         : applyFigureProperties

function [totaldist,vinstant,subdist] = quadspline(x,fx,tinstant,t1,t2)

% Initialize variables
N = numel(x) - 1;
b = zeros(3*N,1);   % column matrix of knowns
A = zeros(3*N,3*N); % coefficient matrix

%% Generate Coefficients of Splines

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

%% Plot

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
plot(x,fx,'ro')

% Axis properties
set(gca,'TickLabelInterpreter','latex')
title('Quadratic Spline Fit')
xlabel('Time ($s$)');
ylabel('Distance ($m$)')

%% Total Distance Traveled

% Initialize variables
totaldist = 0;
k = 1;

% Integration of the velocity function over the current segment

for j = 2:N+1 % N+1 data points where N is the number of splines
    I = (c(k)*x(j)^3/3 + c(k+1)*x(j)^2/2 + c(k+2)*x(j)) ...
        - (c(k)*x(j-1)^3/3 + c(k+1)*x(j-1)^2/2 + c(k+2)*x(j-1));
    totaldist = totaldist + I;

    k = k + 3; % increment to next set of coefficients
end

%% <Insert code that calculates the velocity at tinstant.>

% Find the relevant spline that contains tinstant
spline_idx = findSpline(x,tinstant);

% Starting coefficient index
k = 3*spline_idx - 2;

% Calculate the velocity at tinstant using the derivative of the spline
vinstant = c(k)*tinstant^2 + c(k+1)*tinstant + c(k+2);

%% <Insert code that calculates the distance traveled from t1 to t2.>

% Find the relevant spline that contains t1 and t2
spline_idx_t1 = findSpline(x,t1);
spline_idx_t2 = findSpline(x,t2);

% Initialize variables for distance calculation
subdist = 0;
k = 3*spline_idx_t1 - 2;  % Starting coefficient index for t1 segment

% Loop through spline segments from t1 to t2
for j = (spline_idx_t1):(spline_idx_t2 + 1)

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

% Initialize segment_index to default value
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