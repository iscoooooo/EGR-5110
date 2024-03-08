%  Code provided by Dr. Nissenson

clc; clear; close all;

x0 = 1.2; y0 = 0;          % initial position of spaceship (normalized units)
vx0 = 0; vy0 = -1.0493571; % initial velocity of spaceship (normalized units)
t0  = 0; tf = 10;          % initial time, final time (normalized units)
fd  = 0;                   % deceleration coefficient
N   = 10000;               % initial guess for number of time steps
dt  = (tf-t0)/N;           % dt = time step (normalized units)
tol = 0.01;                % tolerance used for termination criteria

% I will adjust the parameters above when grading your assignment

[x,y,vx,vy,t] = myrk4(t0,tf,dt,N,x0,y0,vx0,vy0,fd,tol);
% You will need to write an RK4 solver in a separate M-file named myrk4.m
% x and y are the coordinates of the spaceship at each time step
% vx and vy are the x and y components of the speed at each time step

disp(['The final (x,y) position is: (' num2str(x(end)) ',' num2str(y(end)) ')'])