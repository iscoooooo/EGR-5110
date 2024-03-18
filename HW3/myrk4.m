function [x,y,vx,vy,t] = myrk4(t0,tf,dt,N,x0,y0,vx0,vy0,fd,tol)

% set of ODEs
f = @(t,X) myODEs(t,X,fd);

mu = 1/82.45;                  % ratio of the moon to earth mass
d  = 384400;                   % characteristic length [km]
RE = 6371;                     % Earth radius [km]
RM = 1740;                     % Moon radius [km]

max_iter = 20;
cond     = true;
k        = 0;

while cond && k < max_iter
    k=k+1;
    %% Solver

    % initialize variables
    x  = zeros(N,1); x(1)  = x0;
    y  = zeros(N,1); y(1)  = y0;
    vx = zeros(N,1); vx(1) = vx0;
    vy = zeros(N,1); vy(1) = vy0;
    t  = zeros(N,1);

    % Begin loop
    for i = 1:N-1

        % 1st slope estimates
        k1 = f(t(i), [x(i);vx(i);y(i);vy(i)]);

        % 2nd slope estimates
        k2 = f(t(i) + 0.5*dt, [x(i) + 0.5*k1(1)*dt; vx(i) + 0.5*k1(2)*dt; ...
            y(i) + 0.5*k1(3)*dt; vy(i) + 0.5*k1(4)*dt]);

        % 3rd slope estimates
        k3 = f(t(i) + 0.5*dt, [x(i) + 0.5*k2(1)*dt; vx(i) + 0.5*k2(2)*dt; ...
            y(i) + 0.5*k2(3)*dt; vy(i) + 0.5*k2(4)*dt]);

        % 4th slope estimates
        k4 = f(t(i) + 1*dt, [x(i) + 1*k3(1)*dt; vx(i) + 1*k3(2)*dt; ...
            y(i) + 1*k3(3)*dt; vy(i) + 1*k3(4)*dt]);

        % Update variables
        x(i+1)  = x(i)  + (dt/6)*(1*k1(1) + 2*k2(1) + 2*k3(1) + 1*k4(1));
        vx(i+1) = vx(i) + (dt/6)*(1*k1(2) + 2*k2(2) + 2*k3(2) + 1*k4(2));
        y(i+1)  = y(i)  + (dt/6)*(1*k1(3) + 2*k2(3) + 2*k3(3) + 1*k4(3));
        vy(i+1) = vy(i) + (dt/6)*(1*k1(4) + 2*k2(4) + 2*k3(4) + 1*k4(4));

        t(i+1) = t(i) + dt;
    end

    %% Termination Criteria

    r1_end(k) = ((x(end) + mu)^2 + y(end)^2)^0.5;

    if k ~= 1
        error = abs((r1_end(k) - r1_end(k-1)) / r1_end(k));
        if error < tol
            cond = false;
        end
    end

    N  = 2*N;
    dt = (tf-t0)/N;

end

%% Plotting

r1 = sqrt((x + mu).^2 + y.^2); % radial distance from earth
v  = sqrt(vx.^2 + vy.^2);      % speed

% Figure 1 properties
figure(1)
set(gcf, ...
    'Units','normalized', ...
    'Position',[0.2,0.2,0.5,0.6], ...
    'DefaultTextInterpreter','latex', ...
    'DefaultLegendInterpreter','latex', ...
    'DefaultAxesFontSize',20, ...
    'Color','k', ...
    'InvertHardCopy', 'off')

% Trajectory plot
hold on
h1 = plot(x,y,'m');
h2 = circle(-mu,0,RE/d,'c','c');
h3 = circle(1-mu,0,RM/d,'w','w');
hold off
title('Spacecraft Trajectory')
xlabel('$x$'), ylabel('$y$')
legend([h1,h2,h3],{'Trajectory','Earth','Moon'});
xlim([min(x)*1.1, max(x)*1.1]);
ylim([min(y)*1.1, max(y)*1.1]);
grid on, box on, axis equal

% Iterate over all axes and legends in the figure and adjust properties
allAxes = findall(gcf, 'type', 'axes');
legends = findobj(gcf, 'Type', 'Legend');
for ax = allAxes'
    set(ax, 'Color', 'k', 'XColor', 'w', 'YColor', 'w');
    ax.Title.Color = 'w'; % Set title color to white
    set(legends, 'Color', 'k');          % Set legend background color to black
    set(legends, 'TextColor', 'w');      % Set legend text color to white
    set(legends, 'EdgeColor', 'w');      % Set legend box edge color to white
end

% Figure 2 properties
figure(2)
set(gcf, ...
    'Units','normalized', ...
    'Position',[0.3,0.2,0.5,0.6], ...
    'DefaultTextInterpreter','latex', ...
    'DefaultLegendInterpreter','latex', ...
    'DefaultAxesFontSize',20, ...
    'Color','k', ...
    'InvertHardCopy', 'off')

% r1 vs. v
hold on
plot(r1,v,'Color','c')
xline(RE/d,'r--','LineWidth',2)
hold off
title('Phase Space')
xlabel('$r_1$'), ylabel('$v$')
legend('Phase','Earth Radius');
xlim([min(r1)*0.9, max(r1)*1.1]);
ylim([min(v)*0.9, max(v)*1.1]);
axis equal, grid on, box on

% Iterate over all axes and legends in the figure and adjust properties
allAxes = findall(gcf, 'type', 'axes');
legends = findobj(gcf, 'Type', 'Legend');
for ax = allAxes'
    set(ax, 'Color', 'k', 'XColor', 'w', 'YColor', 'w');
    ax.Title.Color = 'w'; % Set title color to white
    set(legends, 'Color', 'k');          % Set legend background color to black
    set(legends, 'TextColor', 'w');      % Set legend text color to white
    set(legends, 'EdgeColor', 'w');      % Set legend box edge color to white
end

%% Animation

% <If you want to get fancy, you can create an animation of the spaceship flying around
% the earth and moon for extra credit (see below).>

end

function dXdt = myODEs(t,X,fd)

% ratio of moon to earth mass
mu = 1/82.45;

% unpack states
x  = X(1);
vx = X(2);
y  = X(3);
vy = X(4);

% distances from earth/moon to spacecraft
r1 = sqrt((x + mu)^2 + y^2);
r2 = sqrt((x-(1-mu))^2 + y^2);

% accelerations
ax = 2*vy + x - (1-mu)*(x+mu)/r1^3 - mu*(x-(1-mu))/r2^3 - fd*vx;
ay = -2*vx + y - (1-mu)*y/r1^3 - mu*y/r2^3 - fd*vy;

% construct derivative of state vector
dXdt = [vx; ax; vy; ay];

end

function p = circle(x, y, radius, faceColor, edgeColor)
% circle Plots a circle with specified center, radius, and color.
%
% Inputs:
%   x         - x-coordinate of the circle's center
%   y         - y-coordinate of the circle's center
%   radius    - Radius of the circle
%   faceColor - Color of the circle's face
%   edgeColor - Color of the circle's edge

% Calculate the position of the rectangle that will be drawn as a circle
% The position is defined as [xLeft yBottom width height]
position = [x-radius, y-radius, 2*radius, 2*radius];

% Create the circle with specified properties
rectangle('Position', position, ...
    'Curvature', [1 1], ...
    'FaceColor', faceColor, ...
    'EdgeColor', edgeColor);

% Create a proxy object for the circle in the legend
p = plot(NaN,NaN,edgeColor, 'Marker', 'o', 'MarkerFaceColor', faceColor, 'LineStyle', 'none');

% % Adjust the axes limits to ensure the circle is fully visible
% axis equal;
% xlim([x-radius*1.1, x+radius*1.1]);
% ylim([y-radius*1.1, y+radius*1.1]);
end