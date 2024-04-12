% Code provided by Dr. Nissenson

clc; clear; close all;

alpha = 10e-6; % thermal diffusivity (m^2/s)
kcond = 10; % thermal conductivity (W / m K)
h = 50; % convection coefficient (W / m^2 K)

dx = 0.001; % node spacing, x-direction (m)
dy = dx; % node spacing, y-direction (m)
Lx = 0.04; % fin length, x-direction (m)
Ly = 0.01; % fin length, y-direction (m)
Lz = 0.20; % fin length, z-direction (m)
Nx = round(Lx/dx) + 1; % number of nodes, x-direction (m)
Ny = round(Ly/dy) + 1; % number of nodes, y-direction (m)
P  = 2*(Ly+Lz); % Fin cros-sectional perimeter
Ac = Ly*Lz;     % Fin cross-sectional area
 
t0 = 0; tf = 10*60; % initial time, final time (seconds)
Nt = 30000; % number of time steps
dt = (tf-t0)/Nt; % dt = time step (seconds)
% Adjust certain parameters to meet the stability criteria.

T = zeros(Nx, Ny, Nt); % Create T array
T(:,:,1) = 40; % Initial condition (Celsius)
Tb = 100; % base temperature (Celsius)
Tinf = 25; % free stream temperature (Celsius)
% Parameters above may be adjusted

lam = (alpha*dt/dx^2); % Fourier number
Bi = h*dx/kcond; % Biot number

disp(['number of nodes, x-direction = ',num2str(Nx)])
disp(['number of nodes, y-direction = ',num2str(Ny)])
disp(['lam = ',num2str(lam)])
disp(['Bi = ',num2str(Bi)])

% Calculate the temperature distribution as a function of time,
% the average temperature at the tip at end of simulation,
% and heat rate into the fin at end of simulation.
[T,Ttipsim,Qfinsim,tss] = calcTvstime(T,Nx,Ny,Nt,lam,kcond,h,dx,dt,Lx,Ly,Lz,Bi,Tb,Tinf);

% Constants
C1 = sqrt(h*P/(kcond*Ac)); C2 = (Tb-Tinf)*sqrt(h*P*kcond*Ac);

% Calculate the temperature at tip using 1D approximation
Ttip1D = Tinf + (Tb - Tinf)*(1/(cosh(C1*Lx) + (h/(C1*kcond)*sinh(C1*Lx))));

disp(['Temp at tip, 1D approximation = ',num2str(Ttip1D)])
disp(['Temp at tip, end of simulation = ', num2str(Ttipsim)])

% Calculate the heat rate into fin using 1D approximation
Qfin1D = C2*((sinh(C1*Lx) + (h/(C1*kcond))*cosh(C1*Lx)) / (cosh(C1*Lx) + (h/(C1*kcond))*sinh(C1*Lx)));

disp(['Heat rate into fin, 1D approximation = ',num2str(Qfin1D)])
disp(['Heat rate into fin, end of simulation = ', num2str(Qfinsim)])

disp(['time to 0.01% steady state = ', num2str(tss/60), ' minutes'])