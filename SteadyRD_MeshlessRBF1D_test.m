function SteadyRD_MeshlessRBF1D_test ()
%-------------------------------------------------------------------------%
% Test of SteadyRD_MeshlessRBF1D code                                     %
%-------------------------------------------------------------------------%
% Author : Hamzah Bakhti -------------- Created : 18 July 2018            %
%-------------------------------------------------------------------------%
% Comparaison of the exact solution and RBF meshless numerical solution   %
% for the following reaction-diffusion equation :                         %
% (E) cos(pi x) d^2u/du^2 + pi sin(pi x) du/dx +                          %
%     exp(x) u = exp(x) sin(pi x) in domain [0,1]                         %
% The boundary conditions are : u(0) = 0 and -u(1) - du(1)/dx = pi        %
% The exact solution is : uex(x) = sin(pi x) in [0,1]                     %
%-------------------------------------------------------------------------%
% This code is distributed under the GNU gpl-3.0 license.                 %
% ------------------------------------------------------------------------%
% Space discretization
xmin = 0; xmax = 1;
nx = 50; dx = (xmax-xmin)/nx;
x = xmin:dx:xmax;
% nx = 50; x = [0,sort(rand(1,nx-1)),1]; % random nodes
% RBF parameter
c = 0.1;
% Source term
f = @(x) exp(x).*sin(pi*x);
% System coefficients
a1 = @(x) cos(pi*x); a2 = @(x) pi*sin(pi*x); a3 = @(x) exp(x);
% Boundary coefficients
beta = [1 -1]; lambda = [0 -1]; g = [0 pi];
% Solution
[u, lhs, rhs] = SteadyRD_MeshlessRBF1D(nx, x, c, a1, a2, a3, f, ...
    beta, lambda, g);
% Exacte solution
uex = sin(pi*x)';
err = abs(uex-u);
% Plot and comparaison
figure
plot(x,u,x,uex,'*r')
xlabel('x-axis')
ylabel('u')
title('Exact and numerical solutions')
legend('Meshless Solution','Exacte Solution')
grid on
figure
plot(x,err,'-s')
xlabel('x-axis')
ylabel('L_1 Erreur')
title('Erreur of RBF meshless method')
grid on