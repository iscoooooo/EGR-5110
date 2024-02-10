clc; clear;

% HW1-Activity given constants
a = 30*pi/180; % [rad]
b = 45*pi/180; % [rad]
P = 3;         % [kN]

% Define coefficient matrix
A = [sin(b) cos(a)   0     0    0    0 1 0; 
     cos(b) sin(a)   0     0    0    0 0 1;
    -sin(b)   0    sin(b)  0    0    0 0 0; 
    -cos(b)   0   -cos(b) -1    0    0 0 0;
       0      0   -sin(b)  0 -cos(a) 0 0 0;
       0      0    cos(b)  0  sin(a) 1 0 0; 
       0   -cos(a)   0     0  cos(a) 0 0 0; 
       0   -sin(a)   0     1 -sin(a) 0 0 0];

% Define vector of constants
 B = [0;
      0; 
     -P; 
      0; 
      0; 
      0; 
      0; 
      0];

% start timer
tic;

% Note: The vector xGE in GaussElim should be a horizontal vector when output.
[XGE, determ] = GaussElim(A,B);

% stop timer
time = toc;

% The transpose operator converts XGE into a column vector.
XGE = XGE';

% For comparison with MATLAB’s solver
XMATLAB = A\B;

compare = [XGE, XMATLAB]; % Combines arrays for side-by-side comparison
fprintf('\n\t  GE \t  MATLAB \n')
disp(compare)
disp(['Elapsed time is ', num2str(time), ' seconds'])

disp(['Your solver calculated ', num2str(determ), ' for the determinant'])
disp(['MATLABs det function gives ', num2str(det(A))])
fprintf('\n')