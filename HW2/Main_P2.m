% need to determine dy analytically 

clc; clear;
dens = 999;             % fluid density, kg/m^3
visc = 1e-3;            % fluid viscosity, N s/m^2
d = 0.25;               % pipe diameter, m
Q = 0.05;               % volumetric flow rate, m^3/s
rough = 0.05e-3;        % absolute roughness, m
fra = 0.008; frb = 0.1; % bounds for bisection method
fr0 = 0.01;             % initial guess for N-R method
fr1 = 0.01; fr2 = 0.02; % initial guesses for secant method
tol = 0.00001;          % tolerance
% I will adjust the parameters above when grading your assignment, but will select values
% that will not cause your solvers to fail.

Re = 4*abs(Q)*dens/(pi*d*visc);

y=@(fr) 1./sqrt(fr) + 2.0*log10((rough/d)/3.7 + 2.51./(Re.*sqrt(fr)));
% dy=@(fr) % I will insert the correct derivative expression when testing your code
           % When you test your N-R solver, you will need to determine dy
           
if(Re >= 4000)
    [frBis, iterBis, frNR, iterNR, frSec, iterSec] = RootFinder_P1(y, dy, fra, frb, fr0, fr1, fr2, tol);

    fzeroVal = fzero(y, [fra,frb]);

    disp(['Bisection: fr = ', num2str(frBis), ' in ', num2str(iterBis), ' iterations.'])
    disp(['Newton-Raphson: fr = ', num2str(frNR), ' in ', num2str(iterNR), ' iterations.'])
    disp(['Secant: fr = ', num2str(frSec), ' in ', num2str(iterSec), ' iterations.'])
    disp(['fzero: fr = ', num2str(fzeroVal)])
else
    disp('Using this program is not appropriate')
end