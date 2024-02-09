function [X, determ] = GaussElim(A, B)

n = size(A,1); 
a = [A,B];     
p = 0;       

if ~isnumeric(A) || ~isnumeric(B)
    error("Inputs must be numeric arrays, and not %s",class(n))
elseif size(A,2) ~= size(B,1)
    error("Matrices A and B are not compatible.")
elseif nargin < 2
    error("Not enough input arguments.")
elseif nargin > 2
    error("Too many input arguments.")
end

for k = 1:n-1
    [max_mag,idx]  = myMax(abs(a(k:n,k))); 

    if abs(a(k,k)) < max_mag
        temp     = a(k,:);    
        a(k,:)   = a(k+idx-1,:); 
        a(k+idx-1,:) = temp;    
        p = p + 1;             
    end

    for i = k+1:n
        ratio  = a(i,k)/a(k,k);
        a(i,:) = a(i,:) - ratio*a(k,:);
    end
end

U = a(:,1:n); 
C = a(:,n+1); 

X = zeros(n,1);

X(n) = C(n)/U(n,n); 

for i = n-1:-1:1   
    X(i) = (C(i) - U(i,i+1:n)*X(i+1:n)) / U(i,i);
end

max_error = myMax(A*X - B);

tol   = 1e-6; 
if max_error > tol
    fprintf("Round-off error is significant. Solution may be incorrect.");
end

X = X'; 

product = 1;
for i = 1:n
    product = product*U(i,i);
    if i == n
        determ_U = product;
    end
end

determ = (-1)^p*determ_U;

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