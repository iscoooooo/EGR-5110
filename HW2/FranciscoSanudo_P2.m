function Q = FranciscoSanudo_P2(dens, visc, rough, tol)

if ~isnumeric(dens) || ~isnumeric(visc) || ~isnumeric(rough) || ~isnumeric(tol)
    error('All inputs must be numeric.');
end

if ~isscalar(dens) || ~isscalar(visc) || ~isscalar(rough) || ~isscalar(tol)
    error('All inputs must be scalar.');
end

if dens <= 0 || visc <= 0 || rough <= 0 || tol <= 0
    error('Density, viscosity, roughness, and tolerance must be positive.');
end

L = [200;100;360;200;100;200;300;450]; 

D = [123.4;158.6;123.4;123.4;176.2;96.8;123.4;109.8]/1000;

Q0 = [1200;800;300;900;100;300;300;400]/3600; 

delta = 0.01;

parameters = struct('Lengths',L,'Diameters',D,'Density',dens,'Viscosity',visc, ...
    'Roughness',rough);

myEquations = @(Q) myFunc(Q,parameters);

[Q,~,QPath,normPath] = multivariateRootFinder(myEquations,Q0,delta,tol);

Q1_path = QPath(1,:);

plotResults(Q1_path, normPath)

end

function F = myFunc(Q,parameters)

L     = parameters.Lengths;  
D     = parameters.Diameters; 
dens  = parameters.Density;  
visc  = parameters.Viscosity;  
rough = parameters.Roughness;
g     = 9.81;                

Re = 4.*abs(Q)*dens./(pi.*D.*visc);

fr = (-1.8.*log10(((rough./D)./3.7).^1.11 + 6.9./Re)).^-2;

F = [Q(1) + Q(2) - 2000/3600; ...                    

Q(3) + Q(7) + Q(8) - 1000/3600; ...              

-Q(4) - Q(5) - Q(7) + 1300/3600; ...               

-Q(1) + Q(5) + Q(6) + 800/3600; ...               

-Q(6) - Q(8) + 700/3600; ...                     

-(fr(1)*8*L(1)/(pi^2*g*D(1)^5))*Q(1)*abs(Q(1)) + ... 
(fr(2)*8*L(2)/(pi^2*g*D(2)^5))*Q(2)*abs(Q(2)) + ...
(fr(4)*8*L(4)/(pi^2*g*D(4)^5))*Q(4)*abs(Q(4)) - ...
(fr(5)*8*L(5)/(pi^2*g*D(5)^5))*Q(5)*abs(Q(5));  ...

-(fr(3)*8*L(3)/(pi^2*g*D(3)^5))*Q(3)*abs(Q(3)) - ... 
(fr(4)*8*L(4)/(pi^2*g*D(4)^5))*Q(4)*abs(Q(4)) + ...
(fr(7)*8*L(7)/(pi^2*g*D(7)^5))*Q(7)*abs(Q(7));  ...

(fr(5)*8*L(5)/(pi^2*g*D(5)^5))*Q(5)*abs(Q(5)) - ... 
(fr(6)*8*L(6)/(pi^2*g*D(6)^5))*Q(6)*abs(Q(6)) - ...
(fr(7)*8*L(7)/(pi^2*g*D(7)^5))*Q(7)*abs(Q(7)) + ...
(fr(8)*8*L(8)/(pi^2*g*D(8)^5))*Q(8)*abs(Q(8))];

end

function [X,iter,XPath,normPath] = multivariateRootFinder(FUN,X0,delta,tol)

if ~isa(FUN, 'function_handle')
    error('multivariateRootFinder:InvalidInput', 'FUN must be a function handle.');
end

if ~isnumeric(X0) || ~isvector(X0)
    error('multivariateRootFinder:InvalidInput', 'X0 must be a numeric vector.');
end

if ~isnumeric(delta) || numel(delta) ~= 1 || delta == 0
    error('multivariateRootFinder:InvalidInput', 'delta must be a non-zero numeric scalar.');
end

try
    F0 = FUN(X0);
catch ME
    error('multivariateRootFinder:FunctionEvaluationError', 'Error evaluating FUN at X0: %s', ME.message);
end

if ~isnumeric(F0) || ~isvector(F0)
    error('multivariateRootFinder:InvalidFunctionOutput', 'FUN must return a numeric vector.');
end


max_iterations = 1000;  
converged      = false; 
norm           = 1;     
k              = 0;     
X              = X0;   

while k < max_iterations && ~converged

    k = k + 1;

    F = FUN(X(:,k));

    norm(k+1) = sqrt(F'*F);

    if norm(k+1) < tol
        converged = true;
    end

    J = Jacobian(FUN,X(:,k),delta);

    dX = -J\F;

    X(:,k+1) = X(:,k) + dX;

    if any(isnan(X(:,k+1))) || any(isinf(X(:,k+1)))
        error('multivariateRootFinder:NumericalInstability', 'The computation resulted in NaN or Inf.');
    end

end

if ~converged
    warning('multivariateRootFinder:NonConvergence', 'The function did not converge to a solution.');
end

XPath    = X(:,:);  
normPath = norm(:);  
iter     = k-1;       
X        = X(:,end)'; 

end

function J = Jacobian(FUN,X0,delta)

F = FUN;
m = numel(X0);   
n = numel(F(X0)); 

J = zeros(m,n);

for k = 1:n

    X_perturbed = X0;
    X_perturbed(k) = X_perturbed(k) + delta;

    try
        F_perturbed = FUN(X_perturbed);
    catch ME
        error('Jacobian:FunctionEvaluationError', 'Error evaluating FUN at X_perturbed: %s', ME.message);
    end

    if numel(F_perturbed) ~= n
        error('Jacobian:InconsistentOutput', 'FUN must return outputs of consistent size.');
    end

    J(:,k) = (F(X_perturbed) - F(X0)) / delta;

end

end

function plotResults(Q1_values, norms)

figure;

set(gcf,'units','normalized','position', [0, 0, .4, .5], ...
    'DefaultTextInterpreter','Latex');
movegui(gcf,'center')

subplot(2, 1, 1);
hold on
plot(Q1_values);
yline(Q1_values(end),'k--')
hold off
title('\textbf{$Q_1$ vs. Iteration Number}','FontSize',18);
xlabel('Iteration','FontSize',18);
ylabel('$Q_1 \,\mathrm{\frac{m^3}{s}}$','FontSize',18);
ylim([0, max(Q1_values) + .10*max(Q1_values)])
legend('$Q_1$ path','Solution','interpreter','latex')

subplot(2, 1, 2);
plot(norms);
title('\textbf{Euclidean Norm vs. Iteration Number}','FontSize',18);
xlabel('Iteration','FontSize',18);
ylabel('$\|\mathbf{F}\|_e$','FontSize',18);

end