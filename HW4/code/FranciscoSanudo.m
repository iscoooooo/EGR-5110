function [T,Ttipsim,Qfinsim,tss] = calcTvstime(T,Nx,Ny,Nt,lam,kcond,h,dx,dt,Lx,Ly,Lz,Bi,Tb,Tinf)

t = zeros(Nt,1);

for k = 1:Nt
    for i = 2:Nx-1
        for j = 2:Ny-1
            T(i,j,k+1) = lam*(T(i-1,j,k) + T(i,j-1,k) + T(i+1,j,k) + T(i,j+1,k)) + (1-4*lam)*T(i,j,k);
        end
    end

    for j = 2:Ny-1
        T(Nx,j,k+1) = lam*(2*T(Nx-1,j,k) + T(Nx,j+1,k) + T(Nx,j-1,k) + 2*Bi*Tinf) + (1-4*lam-2*Bi*lam)*T(Nx,j,k);
    end

    for j = 1:Ny
        T(1,j,k+1) = Tb;
    end

    for i = 2:Nx-1
        T(i,Ny,k+1) = lam*(2*T(i,Ny-1,k) + T(i+1,Ny,k) + T(i-1,Ny,k) + 2*Bi*Tinf) + (1-4*lam-2*Bi*lam)*T(i,Ny,k);
    end

    for i = 2:Nx-1
        T(i,1,k+1) = lam*(2*T(i,2,k) + T(i+1,1,k) + T(i-1,1,k) + 2*Bi*Tinf) + (1-4*lam-2*Bi*lam)*T(i,1,k);
    end

    T(Nx,Ny,k+1) = 2*lam*(T(Nx-1,Ny,k) + T(Nx,Ny-1,k) + 2*Bi*Tinf) + (1-4*lam-4*Bi*lam)*T(Nx,Ny,k);

    T(Nx,1,k+1) = 2*lam*(T(Nx-1,1,k) + T(Nx,2,k) + 2*Bi*Tinf) + (1-4*lam-4*Bi*lam)*T(Nx,1,k);

    t(k+1) = t(k) + dt;
end

Ttipsim = (1/((Ny-2) + 0.5 + 0.5)).*(0.5*T(Nx,1,end) + 0.5*T(Nx,Ny,end) + sum(T(Nx,2:end-1,end)));

Qfinsim = 0;
for j = 1:Ny
    if j == 1 || j == Ny
        Qfin = kcond*(0.5*dx*Lz)*(T(1,j,end) - T(2,j,end))/dx;
    else
        Qfin = kcond*(dx*Lz)*(T(1,j,end) - T(2,j,end))/dx;
    end
    Qfinsim = Qfinsim + Qfin;
end

converged = false;
tol = 0.01/100;
k = 0;

while ~converged
    k = k + 1;

    Ttipavg = (1/((Ny-2) + 0.5 + 0.5))*(0.5*T(Nx,1,k) + 0.5*T(Nx,Ny,k) + sum(T(Nx,2:end-1,k)));

    error = abs((Ttipsim - Ttipavg)/Ttipsim);

    if error < tol
        converged = true;
    else
        tss = t(k) + dt;
    end
end

T = permute(T, [2, 1, 3]);

response = input('\nDo you want to play the animation? [y/n]: ','s');

if lower(response) == 'y'
    fprintf('\nPlaying animation...\n\n');

    xval = 0:dx:Lx;
    yval = 0:dx:Ly;
    [x,y] = meshgrid(xval,yval);

    f = figure;
    position = [0.2, 0.2, 0.5, 0.6];
    applyFigureProperties(f, position)

    frameskip = 100;

    dT = 5;

    for k = 2:floor((Nt-1)/frameskip):Nt
        [~,h] = contour(x, y, T(:,:,k));

        axis equal;

        h.LevelList = 0:dT:Tb;
        h.ShowText = 'on';

        c = colorbar;
        title(c,'$T$ (${}^{\circ}$C)','interpreter','latex')
        c.TickLabelInterpreter = 'latex';
        colormap('turbo'),

        set(gca,'TickLabelInterpreter','latex')
        xlabel('\textbf{Horizontal Position} ($m$)');
        ylabel('\textbf{Vertical Position} ($m$)')
        tPlot = sprintf('%.2f',t(k)/60);
        title(['Temperature Distribution (${}^{\circ}$C) at $t$ = ', tPlot, ' minutes']);
        axis equal;

        pause(0.01);
    end

elseif lower(response) == 'n'
    fprintf(['\nAnimation skipped. Displaying temperature distribution at' ...
        ' the end of the simulation.\n\n']);

    xval = 0:dx:Lx;
    yval = 0:dx:Ly;
    [x,y] = meshgrid(xval,yval);

    g = figure;
    position = [0.2, 0.2, 0.5, 0.6];
    applyFigureProperties(g, position)

    dT = 5;

    [~,l] = contour(x, y, T(:,:,end));
    axis equal;

    l.LevelList = 0:dT:Tb;
    l.ShowText = 'on';

    cbar = colorbar;
    title(cbar,'$T$ (${}^{\circ}$C)','interpreter','latex')
    cbar.TickLabelInterpreter = 'latex';
    colormap('turbo'),

    axis equal;
    set(gca,'TickLabelInterpreter','latex')
    xlabel('\textbf{Horizontal Position} ($m$)');
    ylabel('\textbf{Vertical Position} ($m$)')
    tPlot = sprintf('%.2f',t(end)/60);
    title(['Temperature Distribution (${}^{\circ}$C) at $t$ = ', tPlot, ' minutes'])
else
    fprintf('\nInvalid input. Animation skipped.\n\n');
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