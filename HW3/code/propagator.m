function [t, X] = propagator(funcODE, X0, N, dt)

% Solution array sizing
numRows = N;
numCols = numel(X0);

% Initialize variables
t = zeros(N,1);
X = zeros(numRows,numCols);

% Set initial condtions
X(1,:) = X0;

% Begin loop
for i = 1:N-1
    % 1st slope estimates
    k1 = funcODE(t(i), X(i,:));

    % 2nd slope estimates
    k2 = funcODE(t(i) + 0.5*dt, X(i,:) + 0.5*k1*dt);

    % 3rd slope estimates
    k3 = funcODE(t(i) + 0.5*dt, X(i,:) + 0.5*k2*dt);

    % 4th slope estimates
    k4 = funcODE(t(i) + 1*dt, X(i,:) + 1*k3*dt);

    % Update state vector
    X(i+1,:) = X(i,:) + (dt/6)*(1*k1 + 2*k2 + 2*k3 + 1*k4);

    % Update time
    t(i+1) = t(i) + dt;
end

end