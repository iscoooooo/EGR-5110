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
% - T       : Temperature distribution array [Nx x Ny x Nt] (Celsius)
% - Tipsim  : Average temperature at the tip at end of simulation (Celsius)
% - Qinfsim : Heat rate into fin at end of simulation (W)
% - tss     : Time to steady-state temperature (seconds)
%
% OTHER
% .m files required              : MAIN.m (calling script)
% Files required (not .m)        : none
% Built-in MATLAB functions used : sum, zeros, abs, permute, meshgrid, floor
% User-defined functions         : applyFigureProperties

function [T,Ttipsim,Qfinsim,tss] = calcTvstime(T,Nx,Ny,Nt,lam,kcond,h,dx,dt,Lx,Ly,Lz,Bi,Tb,Tinf)
%% Compute Temperature Distribution

% Initialize time array
t = zeros(Nt,1);

% Calculate the temperature distribution at each time step
for k = 1:Nt
    % Interior Nodes
    for i = 2:Nx-1
        for j = 2:Ny-1
            T(i,j,k+1) = lam*(T(i-1,j,k) + T(i,j-1,k) + T(i+1,j,k) + T(i,j+1,k)) + (1-4*lam)*T(i,j,k);
        end
    end

    % Right Boundary Nodes
    for j = 2:Ny-1
        T(Nx,j,k+1) = lam*(2*T(Nx-1,j,k) + T(Nx,j+1,k) + T(Nx,j-1,k) + 2*Bi*Tinf) + (1-4*lam-2*Bi*lam)*T(Nx,j,k);
    end

    % Left Boundary Nodes
    for j = 1:Ny
        T(1,j,k+1) = Tb;
    end

    % Top Boundary Nodes
    for i = 2:Nx-1
        T(i,Ny,k+1) = lam*(2*T(i,Ny-1,k) + T(i+1,Ny,k) + T(i-1,Ny,k) + 2*Bi*Tinf) + (1-4*lam-2*Bi*lam)*T(i,Ny,k);
    end

    % Lower Boundary Nodes
    for i = 2:Nx-1
        T(i,1,k+1) = lam*(2*T(i,2,k) + T(i+1,1,k) + T(i-1,1,k) + 2*Bi*Tinf) + (1-4*lam-2*Bi*lam)*T(i,1,k);
    end

    % Top Right Corner
    T(Nx,Ny,k+1) = 2*lam*(T(Nx-1,Ny,k) + T(Nx,Ny-1,k) + 2*Bi*Tinf) + (1-4*lam-4*Bi*lam)*T(Nx,Ny,k);

    % Bottom Right Corner
    T(Nx,1,k+1) = 2*lam*(T(Nx-1,1,k) + T(Nx,2,k) + 2*Bi*Tinf) + (1-4*lam-4*Bi*lam)*T(Nx,1,k);

    % Update time
    t(k+1) = t(k) + dt;
end

%% Average Temperature @ Tip, Heat Rate, & Time to Steady-State

% Calculate the average temperature at the tip at the last time step
Ttipsim = (1/((Ny-2) + 0.5 + 0.5)).*(0.5*T(Nx,1,end) + 0.5*T(Nx,Ny,end) + sum(T(Nx,2:end-1,end)));

% Calculate the heat rate into the fin at the last time step
Qfinsim = 0;
for j = 1:Ny
    if j == 1 || j == Ny
        Qfin = kcond*(0.5*dx*Lz)*(T(1,j,end) - T(2,j,end))/dx;
    else
        Qfin = kcond*(dx*Lz)*(T(1,j,end) - T(2,j,end))/dx;
    end
    Qfinsim = Qfinsim + Qfin;
end

% Calculate the time needed to reach 0.01% steady state at the tip.
converged = false;
tol = 0.01/100;
k = 0;

while ~converged
    k = k + 1;

    % Average temperature at the tip (at the current time step)
    Ttipavg = (1/((Ny-2) + 0.5 + 0.5))*(0.5*T(Nx,1,k) + 0.5*T(Nx,Ny,k) + sum(T(Nx,2:end-1,k)));

    % Compute relative error
    error = abs((Ttipsim - Ttipavg)/Ttipsim);

    % Check if error is less than the specified tolerance
    if error < tol
        converged = true;
    else
        tss = t(k) + dt;
    end
end

%% Animation

% Transpose Temperature array across all time steps (swap x and y dim)
T = permute(T, [2, 1, 3]);

% Ask user if they want to play the animation
response = input('\nDo you want to play the animation? [y/n]: ','s');

% Check response
if lower(response) == 'y'
    % User wants to play the animation
    fprintf('\nPlaying animation...\n\n');

    % Create 2-D grid coodinates
    xval = 0:dx:Lx;
    yval = 0:dx:Ly;
    [x,y] = meshgrid(xval,yval);

    % Create figure and apply figure properties
    f = figure;
    position = [0.2, 0.2, 0.5, 0.6];
    applyFigureProperties(f, position)

    % Frames to skip
    frameskip = 100;

    % Increment for temperature contours
    dT = 5;

    % Animation loop
    for k = 2:floor((Nt-1)/frameskip):Nt
        % Plot contour at current time-step
        [~,h] = contour(x, y, T(:,:,k));

        axis equal;

        % Contour properties
        h.LevelList = 0:dT:Tb;
        h.ShowText = 'on';

        % Colorbar properties
        c = colorbar;
        title(c,'$T$ (${}^{\circ}$C)','interpreter','latex')
        c.TickLabelInterpreter = 'latex';
        colormap('turbo'),

        % Axis properties
        set(gca,'TickLabelInterpreter','latex')
        xlabel('\textbf{Horizontal Position} ($m$)');
        ylabel('\textbf{Vertical Position} ($m$)')
        tPlot = sprintf('%.2f',t(k)/60);
        title(['Temperature Distribution (${}^{\circ}$C) at $t$ = ', tPlot, ' minutes']);
        axis equal;

        % Adjust animation speed
        pause(0.01);
    end

elseif lower(response) == 'n'
    % User does not want to play the animation
    fprintf(['\nAnimation skipped. Displaying temperature distribution at' ...
        ' the end of the simulation.\n\n']);

    %...Show temperature distribution at the end of the simulation

    % Create 2-D grid coordinates
    xval = 0:dx:Lx;
    yval = 0:dx:Ly;
    [x,y] = meshgrid(xval,yval);

    % Create figure and apply figure properties
    g = figure;
    position = [0.2, 0.2, 0.5, 0.6];
    applyFigureProperties(g, position)

    % Increment for temperature contours
    dT = 5;

    % Plot contours
    [~,l] = contour(x, y, T(:,:,end));
    axis equal;

    % Contour properties
    l.LevelList = 0:dT:Tb;
    l.ShowText = 'on';

    % Colorbar properties
    cbar = colorbar;
    title(cbar,'$T$ (${}^{\circ}$C)','interpreter','latex')
    cbar.TickLabelInterpreter = 'latex';
    colormap('turbo'),

    % Axis properties
    axis equal;
    set(gca,'TickLabelInterpreter','latex')
    xlabel('\textbf{Horizontal Position} ($m$)');
    ylabel('\textbf{Vertical Position} ($m$)')
    tPlot = sprintf('%.2f',t(end)/60);
    title(['Temperature Distribution (${}^{\circ}$C) at $t$ = ', tPlot, ' minutes'])
else
    % Invalid input
    fprintf('\nInvalid input. Animation skipped.\n\n');
end

end

%-------------------------------------------------------------------------

function applyFigureProperties(figHandle, position)
set(figHandle, ...
    'Units', 'normalized', ...
    'Position', position, ...
    'DefaultTextInterpreter', 'latex', ...
    'DefaultLegendInterpreter', 'latex', ...
    'DefaultAxesFontSize', 14);
end