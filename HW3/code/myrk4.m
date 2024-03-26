function [x,y,vx,vy,t] = myrk4(t0,tf,dt,N,x0,y0,vx0,vy0,fd,tol)

% Physical parameters
mu    = 1/82.45;     % Ratio of the moon to earth mass
G     = 6.674e-11;   % [Nm^2kg^-2] Gravitational constant
RE    = 6371;        % [km]        Earth radius
RM    = 1740;        % [km]        Moon radius
massE = 5.972e24;    % [kg]        Earth mass
massM = 7.348e22;    % [kg]        Moon mass
d     = 384400;      % [km]        Characteristic length
ct    = 375203.2339; % [sec]       Characteristic time
cv    = d/ct;        % [km/s]      Characteristic speed

% Set of ODEs
f = @(t,X) myODEs(t,X,mu,fd);

% Initial condtions
X0 = [x0, y0, vx0, vy0];

% Loop conditions
max_iter   = 20;
cond       = true;
k          = 0;
r1_end_old = NaN;

% Begin loop
while cond && k < max_iter

    % Increment counter
    k = k+1;

    %% Solver
    [t,X] = propagator(f, X0, N, dt);

    % Extract position & velocity components
    x = X(:,1); vx = X(:,3);
    y = X(:,2); vy = X(:,4);

    % % check for collision with earth and moon (could place in propagator)
    % r1_min = min(((x + mu).^2 + y.^2).^0.5);
    % r2_min = max(((x - 1 + mu).^2 + y.^2).^0.5);
    %
    % if r1_min < RE
    %     fprintf("Collision with Earth.")
    %     break
    % elseif r2_min < RM
    %     fprintf("Collision with Moon.")
    %     break
    % end

    %% Termination Criteria

    % Get final value of r1
    r1_end = ((x(end) + mu)^2 + y(end)^2)^0.5;

    if ~isnan(r1_end_old) % check if previous r1 value has been set
        error = abs((r1_end - r1_end_old) / r1_end);
        if error < tol
            cond = false;
        end
    end

    % Update previous r1 value for next iteration
    r1_end_old = r1_end;

    % Adjust N and dt for next iteration
    N  = 2*N;
    dt = (tf-t0)/N;

end

%% Plotting

r1 = sqrt((x + mu).^2 + y.^2); % Radial distance from earth
v  = sqrt(vx.^2 + vy.^2);      % Speed

% Figure 1
figure1 = figure(1);
applyFigureProperties(figure1, [0.2, 0.2, 0.5, 0.6]);

% Trajectory plot
hold on;
plot(x, y, 'm','LineWidth',2);
circle(-mu, 0, RE/d, 'c', 'c');
circle(1-mu, 0, RM/d, 'w', 'w');
hold off;

title('Spacecraft Trajectory');
xlabel('$x$'), ylabel('$y$');
legend('Trajectory', 'Earth', 'Moon');
% xlim([min(x)*1.1, max(x)*1.1]);
% ylim([min(y)*1.1, max(y)*1.1]);

applyAxisAndLegendProperties(figure1);

% Figure 2
figure2 = figure(2);
applyFigureProperties(figure2, [0.3, 0.2, 0.5, 0.6]);

% r1 vs. v plot
hold on;
plot(r1, v, 'Color', 'c','LineWidth',2);
xline(RE/d, 'r--', 'LineWidth', 2);
hold off;

title('Phase Space');
xlabel('$r_1$'), ylabel('$v$');
legend('Phase', 'Earth Radius');

applyAxisAndLegendProperties(figure2);


%% Animation

% Ask user if they want to play the animation
response = input('Do you want to play the animation? [y/n]: ','s');

% Check the user's response
if lower(response) == 'y'
    % User wants to play the animation
    fprintf('\nPlaying animation...\n\n');

    % Create figure and apply figure properties
    figure3 = figure(3);
    applyFigureProperties(figure3, [0.2, 0.2, 0.5, 0.6]);

    hold on;
    % Static parts: Earth and Moon
    circle(-mu, 0, RE/d, 'c', 'c');
    circle(1-mu, 0, RM/d, 'w', 'w');

    % Spacecraft spacecraft & trajectory initialization
    hSpacecraft = plot(nan, nan, 'mo', 'MarkerFaceColor', 'm');
    trajectoryPlot = plot(x(1), y(1), 'm','HandleVisibility','off');

    % Set plot limits and labels
    % xlim([min(x)*1.1, max(x)*1.1]);
    % ylim([min(y)*1.1, max(y)*1.1]);
    xlabel('$x$'), ylabel('$y$');
    title('Spacecraft Trajectory Animation');
    legend('Earth', 'Moon','Spacecraft');

    % Apply axis and legend properties after all plot commands
    applyAxisAndLegendProperties(figure3);

    % Pre-rendering the animation with reduced frames
    numPoints = length(x);
    frameSkip = 20; % Adjust frameSkip accordingly
    frames(ceil(numPoints/frameSkip)) = struct('cdata',[],'colormap',[]);

    j = 1;
    for k = 1:frameSkip:numPoints
        % Update trajectory plot to include points up to the current one
        set(trajectoryPlot, 'XData', x(1:k), 'YData', y(1:k));
        % Update spacecraft position
        set(hSpacecraft, 'XData', x(k), 'YData', y(k));
        drawnow;
        frames(j) = getframe(figure3);
        j = j + 1;
    end
    
    hold off;
    
elseif lower(response) == 'n'
    % User does not want to play the animation
    fprintf('\nAnimation skipped.\n\n');
else
    % Invalid input
    fprintf('\nInvalid input. Animation skipped.\n\n');
end

end

%% Sub-functions

%------------------------------------------------------------------------

function dXdt = myODEs(t,X,mu,fd)

% Unpack states
x  = X(1);
y  = X(2);
vx = X(3);
vy = X(4);

% Radial distances from the two primary bodies to spacecraft
r1 = sqrt((x + mu)^2 + y^2);
r2 = sqrt((x-(1-mu))^2 + y^2);

% Accelerations
ax = 2*vy + x - (1-mu)*(x+mu)/r1^3 - mu*(x-(1-mu))/r2^3 - fd*vx;
ay = -2*vx + y - (1-mu)*y/r1^3 - mu*y/r2^3 - fd*vy;

% Construct derivative of state vector
dXdt = [vx, vy, ax, ay];

end

%------------------------------------------------------------------------

function [t, X] = propagator(funcODE, X0, N, dt)

% Solution array sizing
numRows = N;
numCols = numel(X0);

% Initialize variables
t = zeros(N,1);
X = zeros(numRows,numCols);

% Set initial condtions
X(1,:) = X0;

% Begin loop
for i = 1:N-1
    % 1st slope estimates
    k1 = funcODE(t(i), X(i,:));

    % 2nd slope estimates
    k2 = funcODE(t(i) + 0.5*dt, X(i,:) + 0.5*k1*dt);

    % 3rd slope estimates
    k3 = funcODE(t(i) + 0.5*dt, X(i,:) + 0.5*k2*dt);

    % 4th slope estimates
    k4 = funcODE(t(i) + 1*dt, X(i,:) + 1*k3*dt);

    % Update state vector
    X(i+1,:) = X(i,:) + (dt/6)*(1*k1 + 2*k2 + 2*k3 + 1*k4);

    % Update time
    t(i+1) = t(i) + dt;
end

end

%------------------------------------------------------------------------

function applyFigureProperties(figHandle, position)
set(figHandle, ...
    'Units', 'normalized', ...
    'Position', position, ...
    'DefaultTextInterpreter', 'latex', ...
    'DefaultLegendInterpreter', 'latex', ...
    'DefaultAxesFontSize', 20, ...
    'Color', 'k', ...
    'InvertHardCopy', 'off');
end

%------------------------------------------------------------------------

function applyAxisAndLegendProperties(figHandle)
% Iterate over all axes in the figure and adjust properties
allAxes = findall(figHandle, 'type', 'axes');
for ax = allAxes'
    set(ax, ...
        'Color', 'k', ...
        'XColor', 'w', ...
        'YColor', 'w', ...
        'GridColor','w', ...
        'MinorGridColor','w');
    ax.Title.Color = 'w'; % Set title color to white
end

% Apply axes-specific commands
grid(ax, 'on'), box(ax, 'on'), axis(ax, 'equal');

% Adjust legend properties
legends = findobj(figHandle, 'Type', 'Legend');
for lg = legends'
    set(lg, 'Color', 'k', 'TextColor', 'w', 'EdgeColor', 'w');
end

end

%------------------------------------------------------------------------

function p = circle(x, y, radius, faceColor, edgeColor)
% Plots a circle with specified center, radius, and color.

% Calculate the position of the rectangle that will be drawn as a circle
% The position is defined as [xLeft yBottom width height]
position = [x-radius, y-radius, 2*radius, 2*radius];

% Create the circle with specified properties
rectangle('Position', position, ...
    'Curvature', [1 1], ...
    'FaceColor', faceColor, ...
    'EdgeColor', edgeColor);

% Create a proxy object for the circle in the legend
p = plot(NaN,NaN,edgeColor, 'Marker', 'o', 'MarkerFaceColor', faceColor, ...
    'LineStyle', 'none');

end