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

%% Plot Lagrange Point Trajectories (L1-L5)

d     = 384400;      % [km] Characteristic length
RE    = 6371;        % [km] Earth radius
RM    = 1740;        % [km] Moon Radius
fd = 0;              % Deceleration coefficient

t0  = 0; tf = 34;  % intial and final time
N   = 64000;       % number of steps
dt  = (tf-t0)/N;   % time step

% Number of trajectories
numTrajectories = 5;

% Lagrange point initial states:
L_x0 = [0.837023544523, 1.155597402589, -1.005053470159, 0.5-mu, ...
    0.5-mu];
L_y0 = [0, 0, 0, sqrt(3)/2, -sqrt(3)/2];
L_vx0 = [0, 0, -.001, 0, 0];
L_vy0 = [0, 0, -.108, 0.01, 0.01];

% Label & Color
labels = ['L1','L2','L3','L4','L5'];
colors = ['r','g','b','m','w'];

% Initialize cell arrays for outputs
Lx = cell(1, numTrajectories);
Ly = cell(1, numTrajectories);
Lvx = cell(1, numTrajectories);
Lvy = cell(1, numTrajectories);
Lt = cell(1, numTrajectories);

for i = 1:numTrajectories
    % Extract the i-th set of initial conditions
    x0 = L_x0(i);
    y0 = L_y0(i);
    vx0 = L_vx0(i);
    vy0 = L_vy0(i);

    % ODE Options
    options = odeset('RelTol',1e-12,'AbsTol',1e-12);

    % Solve for the i-th trajectory
    [t,X] = ode45(@(t,X) myODEs(t,X,mu,fd), [t0,tf], [x0, y0, vx0, vy0]', options);

    % Store the results
    Lx{i}  = X(:,1);
    Ly{i}  = X(:,2);
    Lvx{i} = X(:,3);
    Lvy{i} = X(:,4);
    Lt{i} = t;
end

% Create figure for trajectory plots
figure4 = figure(4);
applyFigureProperties(figure4, [0.2, 0.2, 0.5, 0.6]);

hold on;

% Plot Earth and Moon - Static parts
circle(-mu, 0, RE/d, 'c', 'c');  % Earth
circle(1-mu, 0, RM/d, 'w', 'w');  % Moon

% Loop over trajectories and plot each
for i = 1:numTrajectories
    plot(Lx{i}, Ly{i}, 'LineWidth',2,'Color',colors(i),'DisplayName',labels(i));
end

% Enhance the plot
title('Lagrange Points');
xlabel('$x$'), ylabel('$y$');
legend('Earth', 'Moon', 'L1', 'L2','L3','L4','L5')
axis equal, grid on
hold off;

applyAxisAndLegendProperties(figure4);