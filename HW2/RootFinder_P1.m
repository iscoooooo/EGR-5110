function [frBis, iterBis, frNR, iterNR, frSec, iterSec] = ...
   RootFinder_P1(y, dy, fra, frb, fr0, fr1, fr2, tol)

%% Bisection Method
i    = 1;   % initialize counter
err  = 1;   % initialize error

% Define search interval
xa   = fra; % lower end
xb   = frb; % upper end

% Begin loop
while err > tol
    xmid(i) = (xa + xb)/2;

    if y(xa)*y(xmid(i)) > 0
        xa = xmid(i);
    else
        xb = xmid(i);
    end
    
    xmid(i+1) = (xa + xb)/2;

    err = abs((xmid(i+1) - xmid(i)) / xmid(i+1));
    i = i + 1;
end

frBis = xmid(end);
iterBis = i;

%% Newton-Raphson Method
i    = 1;   % initialize counter
err  = 1;   % initialize error
x(i) = fr0; % set initial guess

% begin loop
while err > tol
    x(i+1) = x(i) - y(x(i))/dy(x(i));
    err    = abs((x(i+1) - x(i)) / x(i+1));
    i      = i + 1;
end

frNR   = x(end);
iterNR = i-1;

%% Secant Method
clear x

i    = 2;   % initialize counter
err  = 1;   % initialize error

% set initial guesses
x(i-1) = fr1;
x(i)   = fr2;

% begin loop
while err > tol
    x(i+1) = x(i) - y(x(i))*(x(i)-x(i-1)) / (y(x(i))-y(x(i-1)));
    err    = abs((x(i+1) - x(i)) / x(i+1));
    i      = i + 1;
end

frSec   = x(end);
iterSec = i-2;

end