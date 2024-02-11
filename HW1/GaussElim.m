% -----------------------MATLAB Script Information----------------------x
%{
Written by: Francisco Sanudo
Date: 2/1/23
Updated: 2/9/24

PURPOSE
This function solves the linear system AX = B using Gaussian elimination
with partial pivoting. The function also returns the determinant from the
augmented matrix and number of pivots.

REFERENCES
Solving Sets of Linear Algebraic Equations (notes), P. Nissenson

INPUTS
- A : Matrix of size n x n containing n^2 known coefficients
- B : Column matrix of size n x 1 contaning n unknown constants

OUTPUTS
- X      : Solution to AX=B of size n x 1
- determ : Determinant of A

OTHER
.m files required              : none
Files required (not .m)        : none
Built-in MATLAB functions used : size, numel, zeros, abs
User-defined functions         : myMax (nested)
%}

function [X, determ] = GaussElim(A, B)

%% 0.0 Pre-processsing & Error Handling

n = size(A,1); % number of equations/unknowns
a = [A,B];     % augmented matrix
p = 0;         % Set pivot counter

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
    %   It's important to a keep a pivot count any time that a pivot is
    %   made. It will  be important for computing the determinant.

    [max_mag,idx]  = myMax(abs(a(k:n,k))); % return max value and index

    if abs(a(k,k)) < max_mag
        a([k,k+idx-1],:) = a([k+idx-1,k],:); % swap row k with highest mag. row
        p = p + 1;                           % increment pivot counter
    end

    % Forward Elimination pass on all rows i > k:
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
    X(i) = (C(i) - U(i,i+1:n)*X(i+1:n)) / U(i,i);
end

%% 1.2 Round-off Error

% Determine max error between LHS & RHS of the linear system
max_error = myMax(A*X - B);

% Error handling
tol   = 1e-6; % error tolerance
if max_error > tol
    fprintf("Round-off error is significant. Solution may be incorrect.");
end

X = X'; % output as row vector

%% 2.0 Compute Determinant 

% Determinant of a triangular matrix is the product of all its diagonal
% elements. Computing det(U) is as follows:
product = 1;
for i = 1:n
    product = product*U(i,i);
    if i == n
        determ_U = product;
    end
end

% Accounting for the pivots in the forward elimination step, we can find
% the determinant of A:
determ = (-1)^p*determ_U;

%% Custom Max Function

%{
This nested function uses a for loop to iterate between the elements of
the input array 'A' and compares each element to the current maximum
value. After looping through all the elements, the largest number in the
array is returned as 'maximum'. The function also returns the index of
the largest element in the array.
%}

    function [maximum, idx] = myMax(A)
        maximum = A(1);
        idx = 1;
        for j = 2:numel(A)
            if A(j) > maximum
                maximum = A(j);
                idx = j;
            end
        end
    end

end