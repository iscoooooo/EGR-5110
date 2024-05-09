function [xypos,numsteps,numfneval] = calcMaxStudent(f,xi,yi,tol,sigma,beta)

numsteps  = 1;              
numfneval = 0;              
maxiter   = 10000; 
converged = false;

X = [xi; yi];   
xypos = X'; 

while ~converged && numsteps < maxiter
    g = grad(f, X);
    numfneval = numfneval + 4;

    err = euclideanNorm(g);

    if err < tol
        converged = true;
    else
        h = 1; 

        xtemp   = X(1) + g(1)*h;
        ytemp   = X(2) + g(2)*h;

        fcurrent  = f(X(1),X(2));
        numfneval = numfneval + 1;

        while f(xtemp, ytemp) - fcurrent <  sigma*(g'*g)*h
        
            h = beta * h;

            numfneval = numfneval + 1;

            xtemp = X(1) + g(1)*h;
            ytemp = X(2) + g(2)*h;
        end

        X = X + h*g;

        xypos = [xypos; X'];

        numsteps = numsteps + 1;
    end
end

numfneval = numfneval + 1;

x_coords = xypos(:, 1);
y_coords = xypos(:, 2);

max_x_magnitude = myMax(abs(x_coords));
max_y_magnitude = myMax(abs(y_coords));

x_range = [-1, 1] * (1.2 * max_x_magnitude);
y_range = [-1, 1] * (1.2 * max_y_magnitude);

dx = abs((x_range(end) - x_range(1))) / 100;
dy = abs((y_range(end) - y_range(1))) / 100;  

x = x_range(1):dx:x_range(2);
y = y_range(1):dy:y_range(2);

z  = zeros(numel(y),numel(x));

for i = 1:numel(x)
    for j = 1:numel(y)
        z(i,j) = f(x(j), y(i));
    end
end

f = figure;
position = [0.2, 0.2, 0.5, 0.6];
applyFigureProperties(f, position)

hold on; 
ax = gca; 

contour(x,y,z,'ShowText','on');

plot(x_coords,y_coords,'k')
plot(x_coords,y_coords,'ko','MarkerSize',3,'MarkerFaceColor','k')

plot(xi,yi,'bo','Markersize',6,'MarkerFaceColor','b')

plot(x_coords(end),y_coords(end),'mo','Markersize',6,'MarkerFaceColor','m')

set(ax,'TickLabelInterpreter','latex')
title('Gradient Ascent w/ Armijo Condition')
xlabel('$x$')
ylabel('$y$')
legend('','Path','','$x_i$, $y_i$','$x^*$, $y^*$')
grid on

colormap(ax,'turbo');
cb = colorbar;
cb.TickLabelInterpreter = 'latex';
title(cb,'$\rm Level$','Interpreter','latex')

end

function g = grad(f,X)

n = numel(X);         
g = zeros(n, 1);     
delta = 0.00001;       

for i = 1:n
    X_perturbed = X;

    X_perturbed(i) = X_perturbed(i) + delta;

    f_forward = f(X_perturbed(1),X_perturbed(2));

    X_perturbed(i) = X_perturbed(i) - 2 * delta;

    f_backward = f(X_perturbed(1),X_perturbed(2));

    g(i) = (f_forward - f_backward) / (2 * delta);
end

end

function euclidean_norm_value = euclideanNorm(vec)
squared_elements = vec.^2;

sum_of_squares = 0;

for i = 1:numel(vec)
    squared_element = vec(i)^2;

    sum_of_squares = sum_of_squares + squared_element;
end

euclidean_norm_value = sqrt(sum_of_squares);

end

function [maximum, idx] = myMax(A)
maximum = A(1);
idx = 1;         

for j = 2:numel(A)
    if A(j) > maximum
        maximum = A(j);
        idx = j;
    end
end

end

function applyFigureProperties(figHandle, position)
set(figHandle, ...
    'Units', 'normalized', ...
    'Position', position, ...
    'DefaultTextInterpreter', 'latex', ...
    'DefaultLegendInterpreter', 'latex', ...
    'DefaultAxesFontSize', 14);

end