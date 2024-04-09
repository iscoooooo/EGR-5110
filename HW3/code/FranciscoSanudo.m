function [x,y,vx,vy,t] = myrk4(t0,tf,dt,N,x0,y0,vx0,vy0,fd,tol)

mu    = 1/82.45;    
G     = 6.674e-11;  
RE    = 6371;      
RM    = 1740;       
massE = 5.972e24;   
massM = 7.348e22; 
d     = 384400;    
ct    = 375203.2339; 
cv    = d/ct;        

f = @(t,X) myODEs(t,X,mu,fd);

X0 = [x0, y0, vx0, vy0];

max_iter   = 20;
cond       = true;
k          = 0;
r1_end_old = NaN;

while cond && k < max_iter

    k = k+1;

    [t,X] = propagator(f, X0, N, dt);

    x = X(:,1); vx = X(:,3);
    y = X(:,2); vy = X(:,4);

    r1_end = ((x(end) + mu)^2 + y(end)^2)^0.5;

    if ~isnan(r1_end_old) 
        error = abs((r1_end - r1_end_old) / r1_end);
        if error < tol
            cond = false;
        end
    end

    r1_end_old = r1_end;

    N  = 2*N;
    dt = (tf-t0)/N;

end

r1 = sqrt((x + mu).^2 + y.^2); 
v  = sqrt(vx.^2 + vy.^2);      

figure1 = figure(1);
applyFigureProperties(figure1, [0.2, 0.2, 0.5, 0.6]);

hold on;
plot(x, y, 'm','LineWidth',2);
circle(-mu, 0, RE/d, 'c', 'c');
circle(1-mu, 0, RM/d, 'w', 'w');
hold off;

title('Spacecraft Trajectory');
xlabel('$x$'), ylabel('$y$');
legend('Trajectory', 'Earth', 'Moon');

applyAxisAndLegendProperties(figure1);

figure2 = figure(2);
applyFigureProperties(figure2, [0.3, 0.2, 0.5, 0.6]);

hold on;
plot(r1, v, 'Color', 'c','LineWidth',2);
xline(RE/d, 'r--', 'LineWidth', 2);
hold off;

title('Phase Space');
xlabel('$r_1$'), ylabel('$v$');
legend('Phase', 'Earth Radius');

applyAxisAndLegendProperties(figure2);

response = input('Do you want to play the animation? [y/n]: ','s');

if lower(response) == 'y'
    fprintf('\nPlaying animation...\n\n');

    figure3 = figure(3);
    applyFigureProperties(figure3, [0.2, 0.2, 0.5, 0.6]);

    hold on;
    circle(-mu, 0, RE/d, 'c', 'c');
    circle(1-mu, 0, RM/d, 'w', 'w');

    hSpacecraft = plot(nan, nan, 'mo', 'MarkerFaceColor', 'm');
    trajectoryPlot = plot(x(1), y(1), 'm','HandleVisibility','off');

    xlabel('$x$'), ylabel('$y$');
    title('Spacecraft Trajectory Animation');
    legend('Earth', 'Moon','Spacecraft');

    applyAxisAndLegendProperties(figure3);

    numPoints = length(x);
    frames(ceil(numPoints/frameSkip)) = struct('cdata',[],'colormap',[]);

    j = 1;
    for k = 1:frameSkip:numPoints
        set(trajectoryPlot, 'XData', x(1:k), 'YData', y(1:k));
        set(hSpacecraft, 'XData', x(k), 'YData', y(k));
        drawnow;
        frames(j) = getframe(figure3);
        j = j + 1;
    end
    
    hold off;
    
elseif lower(response) == 'n'
    fprintf('\nAnimation skipped.\n\n');
else
    fprintf('\nInvalid input. Animation skipped.\n\n');
end

fprintf("The closest approach to Earth is %.5f (%.2f km)\n\n",min(r1)-RE/d,min(r1)*d-RE)

end

function dXdt = myODEs(t,X,mu,fd)

x  = X(1);
y  = X(2);
vx = X(3);
vy = X(4);

r1 = sqrt((x + mu)^2 + y^2);
r2 = sqrt((x-(1-mu))^2 + y^2);

ax = 2*vy + x - (1-mu)*(x+mu)/r1^3 - mu*(x-(1-mu))/r2^3 - fd*vx;
ay = -2*vx + y - (1-mu)*y/r1^3 - mu*y/r2^3 - fd*vy;

dXdt = [vx, vy, ax, ay];

end

function [t, X] = propagator(funcODE, X0, N, dt)

numRows = N;
numCols = numel(X0);

t = zeros(N,1);
X = zeros(numRows,numCols);

X(1,:) = X0;

for i = 1:N-1
    k1 = funcODE(t(i), X(i,:));

    k2 = funcODE(t(i) + 0.5*dt, X(i,:) + 0.5*k1*dt);

    k3 = funcODE(t(i) + 0.5*dt, X(i,:) + 0.5*k2*dt);

    k4 = funcODE(t(i) + 1*dt, X(i,:) + 1*k3*dt);

    X(i+1,:) = X(i,:) + (dt/6)*(1*k1 + 2*k2 + 2*k3 + 1*k4);

    t(i+1) = t(i) + dt;
end

end

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


function applyAxisAndLegendProperties(figHandle)

allAxes = findall(figHandle, 'type', 'axes');
for ax = allAxes'
    set(ax, ...
        'Color', 'k', ...
        'XColor', 'w', ...
        'YColor', 'w', ...
        'GridColor','w', ...
        'MinorGridColor','w');
    ax.Title.Color = 'w'; 
end

grid(ax, 'on'), box(ax, 'on'), axis(ax, 'equal');

legends = findobj(figHandle, 'Type', 'Legend');
for lg = legends'
    set(lg, 'Color', 'k', 'TextColor', 'w', 'EdgeColor', 'w');
end

end

function p = circle(x, y, radius, faceColor, edgeColor)

position = [x-radius, y-radius, 2*radius, 2*radius];

rectangle('Position', position, ...
    'Curvature', [1 1], ...
    'FaceColor', faceColor, ...
    'EdgeColor', edgeColor);

p = plot(NaN,NaN,edgeColor, 'Marker', 'o', 'MarkerFaceColor', faceColor, ...
    'LineStyle', 'none');

end