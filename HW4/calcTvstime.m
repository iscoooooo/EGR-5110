% Written by: Francisco Sanudo
% Date: 4/4/24
%
% PURPOSE
% calcTvstime solves the 2D Heat Equation to find the temperature
% distribution across a thin fin using an explicit finite-difference method.
%
% REFERENCES
% Solving Partial Differential Equations (notes), P. Nissenson
%
% INPUTS
% - T     : Initialized Temperature Array (Celsius)
% - Tb    : Base Temperature (Celsius)
% - Tinf  : Free-stream Temperature (Celsius)
% - Lx    : Fin length, x-direction (m)
% - Ly    : Fin length, y-direction (m)
% - Lz    : Fin length, z-direction (m)
% - dx    : Node Spacing, x-direction (m)
% - dt    : Node Spacing, y-direction (m)
% - Nx    : Number of nodes, x-direction (m)
% - Ny    : Number of nodes, y-direction (m)
% - Nt    : Number of time steps
% - kcond : Thermal Conductivity (W / m K)
% - h     : Convection Coefficient (W / m^2 K)
% - Bi    : Biot Number
% - lam   : Fourier number
%
% OUTPUTS
% - T       : Temperature distribution
% - Tipsim  : Average temperature at the tip at end of simulation
% - Qinfsim : Heat rate into fin at end of simulation
% - tss     : Time to steady-state temperature
%
% OTHER
% .m files required              :
% Files required (not .m)        :
% Built-in MATLAB functions used :
% User-defined functions         :

function [T,Ttipsim,Qfinsim,tss] = calcTvstime(T,Nx,Ny,Nt,lam,kcond,h,dx,dt,Lx,Ly,Lz,Bi,Tb,Tinf)

% Initialize time array
t = zeros(Nt,1);

% Calculate the temperature distribution at each time step
for k = 1:Nt-1
    %% Interior Nodes
    for i = 2:Nx-1
        for j = 2:Ny-1
            T(i,j,k+1) = lam*(T(i-1,j,k) + T(i,j-1,k) + T(i,j+1,k) + T(i,j+1,k)) + (1-4*lam)*T(i,j,k);
        end
    end

    %% Right Boundary Nodes
    for j = 2:Ny-1
        T(Nx,j,k+1) = lam*(2*T(Nx-1,j,k) + T(Nx,j+1,k) + T(Nx,j-1,k) + 2*Bi*Tinf) + (1-4*lam-2*Bi*lam)*T(Nx,j,k);
    end

    %% Left Boundary Nodes
    for j = 1:Ny
        T(1,j,k+1) = Tb;
    end

    %% Top Boundary Nodes ( derive )


    %% Lower Boundary Nodes ( derive )


    %% Top Right Corner
    T(Nx,Ny,k+1) = 2*lam*(T(Nx-1,Ny,k) + T(Nx,Ny-1,k) + 2*Bi*Tinf) + (1-4*lam-4*Bi*lam)*T(Nx,Ny,k);

    %% Bottom Right Corner ( derive )


    %% Update time
    t(k+1) = t(k) + dt;

end

% Insert code that calculates the average temperature at the tip at the last time step
Ttipsim = 0;

% Insert code that calculates the heat rate into the fin at the last time step
Qfinsim = 0;

% Insert code to determine the time needed to reach 0.01% steady state at the tip.
tss = 0;

% Insert code that creates an animation of the temperature distribution (only plot 100 time steps)

xval = 0:dx:Lx;
yval = 0:dx:Ly;

[x,y] = meshgrid(xval,yval);

[~,h] = contour(x,y,T(:,:,5000)');
axis equal
h.LevelList = 0:5:200;
colorbar