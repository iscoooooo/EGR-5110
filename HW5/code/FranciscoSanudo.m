function [totaldist,vinstant,subdist] = quadspline(x,fx,tinstant,t1,t2)

N = numel(x) - 1;  
b = zeros(3*N,1);  
A = zeros(3*N,3*N);

j = 1; 
k = 1;

for i = 2:2:2*N
    A(i,k:k+2)   = [x(j)^2, x(j), 1];
    b(i)         = fx(j);
    j            = j + 1;           
    A(i+1,k:k+2) = [x(j)^2, x(j), 1];
    b(i+1)       = fx(j);
    k            = k + 3;           
end

j = 1; 
k = 1;

for i = 2*N+2:3*N
    A(i,k:k+4) = [2*x(j), 1, 0, -2*x(j), -1];
    j          = j + 1;
    k          = k + 3;
end

A(1,1) = 1;

c = A\b;

j = 1; 

f = figure;
position = [0.2, 0.2, 0.5, 0.6];
applyFigureProperties(f, position)

hold on 

for i = 1:N
    spline = @(X) c(j)*X.^2 + c(j+1).*X + c(j+2);
    fplot(spline,[x(i) x(i+1)])                 

    j = j + 3; 
end

h = plot(x,fx,'ro');

set(gca,'TickLabelInterpreter','latex')
title('Quadratic Spline Fit')
xlabel('Time ($s$)');
ylabel('Distance ($m$)')
legend(h, 'Data Points')
grid on

totaldist = 0;
k = 1;         


for j = 2:N+1 
    I = (c(k)*x(j)^3/3 + c(k+1)*x(j)^2/2 + c(k+2)*x(j)) ...
        - (c(k)*x(j-1)^3/3 + c(k+1)*x(j-1)^2/2 + c(k+2)*x(j-1));
    totaldist = totaldist + I;

    k = k + 3; 
end

spline_idx = findSpline(x,tinstant);

k = 3*spline_idx - 2;

vinstant = c(k)*tinstant^2 + c(k+1)*tinstant + c(k+2);

spline_idx_t1 = findSpline(x,t1);
spline_idx_t2 = findSpline(x,t2);

subdist = 0;
k = 3*spline_idx_t1 - 2;  

for j = (spline_idx_t1):(spline_idx_t2)

    if j == spline_idx_t1
        I = (c(k)*x(j+1)^3/3 + c(k+1)*x(j+1)^2/2 + c(k+2)*x(j+1)) ...
            - (c(k)*t1^3/3 + c(k+1)*t1^2/2 + c(k+2)*t1);
        subdist = subdist + I;
    elseif j == spline_idx_t2
        I = (c(k)*t2^3/3 + c(k+1)*t2^2/2 + c(k+2)*t2) ...
            - (c(k)*x(j)^3/3 + c(k+1)*x(j)^2/2 + c(k+2)*x(j));
        subdist = subdist + I;
    else 
        I = (c(k)*x(j+1)^3/3 + c(k+1)*x(j+1)^2/2 + c(k+2)*x(j+1)) ...
            - (c(k)*x(j)^3/3 + c(k+1)*x(j)^2/2 + c(k+2)*x(j));
        subdist = subdist + I;
    end

    k = k + 3;
end

end

function idx = findSpline(x, tinstant)

idx = 0;

for i = 1:numel(x)-1
    if tinstant >= x(i) && tinstant < x(i+1)
        idx = i;
        break; 
    end
end

if isempty(idx)
    error('tinstant is out of the range of input time vector x.');
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