clc; clear;

A = [2 0 0 1 1 4 0 0;
     0 0 1 2 1 0 1 1;
     3 3 3 3 3 3 3 3;
     0 0 0 0 1 2 0 0;
     0 0 0 0 0 0 3 1;
     1 1 1 0 0 0 1 2;
     0 5 2 1 0 0 1 2;
     0 2 0 2 0 2 0 2];

B = [21.50;
     24.70;
     70.35;
     7.50;
     13.75;
     13.45;
     28.90;
     27.50];

% start timer
tic;

% Note: The vector xGE in GaussElim should be a horizontal vector when output.
[XGE, determ] = GaussElim(A, B);

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