% -------------MATLAB Script Information-------------x
%
% Written by: Francisco Sanudo
% Date: 2/1/23
% 
% PURPOSE
% This function solves the linear system AX = B using Gaussian elimination
% with partial pivoting. The function also returns the determinant from the
% augmented matrix and number of pivots.  
% 
% REFERENCES
% Solving Sets of Linear Algebraic Equations (notes), P. Nissenson
% 
% INPUTS
% - A : Matrix of size n x n containing n^2 known coefficients
% - B : Column matrix of size n x 1 contaning n unknown constants
%
% OUTPUTS
% - X      : Solution to AX=B of size n x 1
% - determ : Determinant of A
%
% OTHER
% .m files required              : none
% Files required (not .m)        : none
% Built-in MATLAB functions used : size, numel, zeros, abs

function [X, determ] = GaussElim(A, B)
%% 0.0 Pre-processsing & Error Handling

n = size(A,1); % number of equations/unknowns
a = [A,B];     % augmented matrix

if ~isnumeric(A) || ~isnumeric(B)
    error("Inputs must be numeric arrays, and not %s",class(n))
elseif size(A,2) ~= size(B,1)
    error("Matrices A and B are not compatible.")
elseif nargin < 2
    error("Not enough input arguments.")
elseif nargin > 2
    error("Too many input arguments.")
end 

%% 1.0 Gaussian Elimination w/ Partial Pivoting

for k = 1:n-1
    % Partial Pivoting
    %   Search all rows > k for the highest magnitude. If needed, swap row
    %   k with the highest magnitude row.
    [~,idx]  = max(abs(a(k:n,k)));
    temp     = a(k,:);
    a(k,:)   = a(k+idx-1,:);
    a(k+idx-1,:) = temp;

    % Elimination pass on all rows i > k:
    for i = k+1:n
        ratio  = a(i,k)/a(k,k);
        a(i,:) = a(i,:) - ratio*a(k,:);
    end
end

U = a(:,1:n); % upper triangular matrix
C = a(:,n+1); % altered B vector

%% 1.1 Backward Substitution
X = zeros(n,1); % initialize solution array

% Implement backward substitution
X(n) = C(n)/U(n,n); % i=n

for i = n-1:-1:1    % i=n-1,..,1
    total = C(i);
    for j = i+1:n
        total = total - U(i,j)*X(j);
    end
    X(i) = total*(1/U(i,i));
end

X = X'; % output as row vector

%% 1.2 Round-off Error


%% 2.0 Compute Determinant

determ = 0; % placeholder

end