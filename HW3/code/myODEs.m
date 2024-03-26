function dXdt = myODEs(t,X,mu,fd)

% Unpack states
x  = X(1);
y  = X(2);
vx = X(3);
vy = X(4);

% Radial distances from the two primary bodies to spacecraft
r1 = sqrt((x + mu)^2 + y^2);
r2 = sqrt((x-(1-mu))^2 + y^2);

% Accelerations
ax = 2*vy + x - (1-mu)*(x+mu)/r1^3 - mu*(x-(1-mu))/r2^3 - fd*vx;
ay = -2*vx + y - (1-mu)*y/r1^3 - mu*y/r2^3 - fd*vy;

% Construct derivative of state vector
dXdt = [vx, vy, ax, ay]';

end