clc; clear; close all;

%% EDITABLE

customSettings = 'off';      % [on/off] Set to off to choose from different cases

x0 = 0.0295993756503642; y0 = 0;           % initial position of spaceship (normalized units)
vx0 = 0; vy0 = -5.77721960562869;  % initial velocity of spaceship (normalized units)
t0  = 0; tf = 1;            % initial time, final time (normalized units)
fd  = 0;                    % deceleration coefficient
N   = 1000;                 % initial guess for number of time steps
dt  = (tf-t0)/N;            % dt = time step (normalized units)
tol = 0.01;                 % tolerance used for termination criteria
mu = 1/82.45;               % Define mu for 3 body system

%% DO NOT EDIT

if strcmp(customSettings,'off')
    fprintf('Choose from the following cases: \n\n')
    
    fprintf('\t [1]  Infinity shape\n')
    fprintf('\t [2]  Periodic Orbit 1\n')
    fprintf('\t [3]  Periodic Orbit 2\n')
    fprintf('\t [4]  Periodic Orbit 3\n')
    fprintf('\t [5]  Periodic Orbit 4\n')
    fprintf('\t [6]  L1\n')
    fprintf('\t [7]  L2\n')
    fprintf('\t [8]  L3\n')
    fprintf('\t [9]  L4\n')
    fprintf('\t [10] L5\n\n')
    n = input('Selection: ');
    fprintf('\n')

    switch n
        case 1 % Infinity shape
            x0 = 1.2; y0 = 0;
            vx0 = 0; vy0 = -1.0493571;
            t0  = 0; tf = 15;
        case 2 % Periodic Orbit 1
            x0 = 0.879962; y0 = 0;
            vx0 = 0; vy0 = -.38089067106386964470;
            t0  = 0; tf = 19.138746281183026809;
        case 3 % Periodic Orbit 2
            x0 = 0.1003e1; y0 = 0;
            vx0 = 0; vy0 = -0.14465123738451062297e1;
            t0  = 0; tf = 0.12267904265603897140e2;
        case 4 % Periodic Orbit 3
            x0 = 0.994; y0 = 0;
            vx0 = 0; vy0 = -0.21138987966945026683e1;
            t0  = 0; tf = 0.54367954392601899690e1;
        case 5 % Periodic Orbit 4
            x0 = 0.997; y0 = 0;
            vx0 = 0; vy0 = -0.165251217072210773125e1;
            t0  = 0; tf = 0.22929723423442969481e2;
        case 6 % L1
            x0 = 0.837023544523; y0 = 0;
            vx0 = 0; vy0 = 0;
            t0  = 0; tf = 20;
        case 7 % L2
            x0 = 1.155597402589; y0 = 0;
            vx0 = 0; vy0 = 0;
            t0  = 0; tf = 20;
        case 8 % L3
            x0 = -1.005053470159; y0 = 0;
            vx0 = -.001; vy0 = -.108;
            t0  = 0; tf = 20;
        case 9 % L4
            x0 = 0.5-mu; y0 = sqrt(3)/2;
            vx0 = 0; vy0 = 0.01;
            t0  = 0; tf = 20;
        case 10 % L5
            x0 = 0.5-mu; y0 = -sqrt(3)/2;
            vx0 = 0; vy0 = 0.01;
            t0  = 0; tf = 20;
    end
end

[x,y,vx,vy,t] = myrk4(t0,tf,dt,N,x0,y0,vx0,vy0,fd,tol);
% x and y are the coordinates of the spaceship at each time step
% vx and vy are the x and y components of the speed at each time step

fprintf(['The final (x,y) position is: (' num2str(x(end)) ',' num2str(y(end)) ')\n\n'])
