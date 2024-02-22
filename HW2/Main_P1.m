clc; clear;

dens = 999;             % fluid density, kg/m^3
visc = 1e-3;            % fluid viscosity, N s/m^2
d = 0.25;               % pipe diameter, m
Q = 0.05;               % volumetric flow rate, m^3/s
rough = 0.05e-3;        % absolute roughness, m
fra = 0.008; frb = 0.1; % bounds for bisection method
fr0 = 0.01;             % initial guess for N-R method
fr1 = 0.01; fr2 = 0.02; % initial guesses for secant method
tol = 1e-5;             % tolerance
% select values that will not cause your solvers to fail.

Re = 4*abs(Q)*dens/(pi*d*visc); % Reynolds Number

y=@(fr) 1/sqrt(fr) + 2.0*log10((rough/d)/3.7 + 2.51/(Re*sqrt(fr))); % Colebrook Equation
dy=@(fr) - 1/(2*fr^(3/2)) - 251/(100*Re*fr^(3/2)*log(10)*((10*rough)/(37*d) + 251/(100*Re*fr^(1/2)))); % derivative
           
if(Re >= 4000)
    tic
    [frBis, iterBis, frNR, iterNR, frSec, iterSec] = RootFinder_P1(y, dy, fra, frb, fr0, fr1, fr2, tol);
    toc
    
    fzeroVal = fzero(y, [fra,frb]);

    disp(['Bisection: fr = ', num2str(frBis), ' in ', num2str(iterBis), ' iterations.'])
    disp(['Newton-Raphson: fr = ', num2str(frNR), ' in ', num2str(iterNR), ' iterations.'])
    disp(['Secant: fr = ', num2str(frSec), ' in ', num2str(iterSec), ' iterations.'])
    disp(['fzero: fr = ', num2str(fzeroVal)])
else
    disp('Using this program is not appropriate')
end