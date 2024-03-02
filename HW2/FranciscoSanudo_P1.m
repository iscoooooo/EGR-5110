function [frBis, iterBis, frNR, iterNR, frSec, iterSec] = ...
    FranciscoSanudo_P1(y, dy, fra, frb, fr0, fr1, fr2, tol)

if ~isnumeric(fra) || ~isnumeric(frb) || ~isnumeric(fr0) || ...
        ~isnumeric(fr1) || ~isnumeric(fr2) || ~isnumeric(tol)
    error('RootFinder_P1:InvalidInput', 'Arguments must be a numeric value.');
end

if ~isscalar(fra) || ~isscalar(frb) || ~isscalar(fr0) || ...
        ~isscalar(fr1) || ~isscalar(fr2) || ~isscalar(tol)
    error('RootFinder_P1:InvalidInput', 'Arguments must be a scalar.');
end

if isa(y, 'function_handle')
    info = functions(y);
    if ~isfield(info, 'function') && isempty(info.function)
        error('RootFinder_P1:InvalidInput', ['The input y is a function handle,' ...
            ' but not an anonymous function.']);
    end
else
    error('RootFinder_P1:InvalidInput', 'The input y must be an anonymous function.')
end

if isa(dy, 'function_handle')
    info = functions(dy);
    if ~isfield(info, 'function') && isempty(info.function)
        error('RootFinder_P1:InvalidInput', ['The input dy is a function handle, ' ...
            'but not an anonymous function.']);
    end
else
    error('RootFinder_P1:InvalidInput', 'The input dy must be an anonymous function.')
end

[frBis,iterBis] = Bisection(y,fra,frb,tol);

[frNR,iterNR] = newtonRaphson(y,dy,fr0,tol);

[frSec,iterSec] = Secant(y,fr1,fr2,tol);

tol_values  = logspace(-10, -5, 100);  
iter_values = zeros(size(tol_values)); 

for k = 1:numel(tol_values)
    tol_k = tol_values(k);

    [~,i] = Bisection(y,fra,frb,tol_k);

    iter_values(k) = i;
end

figure;

set(gcf,'units','normalized','position', [0, 0, .4, .5], ...
    'DefaultTextInterpreter','Latex');
movegui(gcf,'center')

semilogx(tol_values, iter_values, ...
    'LineStyle', ":","Color","blue",'Marker',".",'MarkerSize',12);

title('\textbf{Bisection Method}: Iterations vs Tolerance','FontSize',15);
xlabel('Tolerance','FontSize',15);
ylabel('Iterations','FontSize',15);
legend('Iterations','fontsize',14,'interpreter','latex')
grid on;

end

function [X,iter] = Bisection(y,a,b,tol)
max_iterations = 1000; 
converged      = false;
i    = 0;               

xa   = a; 
xb   = b; 


while i < max_iterations && ~converged
    i = i + 1;       

    xmid(i) = (xa + xb)/2;

    if y(xa)*y(xmid(i)) > 0
        xa = xmid(i);
    else
        xb = xmid(i);
    end

    xmid(i+1) = (xa + xb)/2;

    err = abs((xmid(i+1) - xmid(i)) / xmid(i+1));
    if err < tol
        converged = true;
    end

    if any(isnan(xmid(i+1))) || any(isinf(xmid(i+1)))
        error('Bisection:NumericalInstability', ['The computation resulted ' ...
            'in NaN or Inf.']);
    end
end

if ~converged
    warning('Bisection:NonConvergence', ['The function did not ' ...
        'converge to a solution.']);
end

X    = xmid(end);
iter = i;        

end

function [X,iter] = newtonRaphson(y,dy,X0,tol)
max_iterations = 1000; 
converged      = false; 
x    = X0;            
i    = 0;               

while i < max_iterations && ~converged
    i = i + 1;

    x(i+1) = x(i) - y(x(i))/dy(x(i));
    
    err    = abs((x(i+1) - x(i)) / x(i+1)); 
    if err < tol
        converged = true;
    end

    if any(isnan(x(i+1))) || any(isinf(x(i+1)))
        error('newtonRaphson:NumericalInstability', ['The computation resulted ' ...
            'in NaN or Inf.']);
    end
end

if ~converged
    warning('newtonRaphson:NonConvergence', ['The function did not ' ...
        'converge to a solution.']);
end

X    = x(end); 
iter = i;     

end

function [X,iter] = Secant(y,X1,X2,tol)
max_iterations = 1000; 
converged      = false; 
i = 1;                 

x(i) = X1;
x(i+1) = X2;

while i < max_iterations && ~converged
    i = i + 1; 

    x(i+1) = x(i) - y(x(i))*(x(i)-x(i-1)) ... 
        / (y(x(i))-y(x(i-1)));

    err = abs((x(i+1) - x(i)) / x(i+1));
    if err < tol
        converged = true;
    end

    if any(isnan(x(i+1))) || any(isinf(x(i+1)))
        error('Secant:NumericalInstability', ['The computation resulted ' ...
            'in NaN or Inf.']);
    end

end

if ~converged
    warning('Secant:NonConvergence', ['The function did not ' ...
        'converge to a solution.']);
end

X    = x(end);
iter = i-1;   

end