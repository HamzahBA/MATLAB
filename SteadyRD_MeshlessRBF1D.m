function [ u, lhs, rhs ] = SteadyRD_MeshlessRBF1D ( nx, x, c, a1, a2, ...
    a3, f, beta, lambda, g )
%-------------------------------------------------------------------------%
% 1D RBF meshless code for steady reaction-diffusion equation             %
%-------------------------------------------------------------------------%
% Author : Hamzah Bakhti ------------------------- Created : 18 July 2018 %
%-------------------------------------------------------------------------%
% Mathematical Model : Steady Reaction-Diffusion Equation                 %
% (E) a1(x) d^2u/du^2 + a2(x) du/dx + a3(x) u = f(x) in [xmin,xmax]       %
% beta_1 u + lambda_1 du/dx = g_1 for x = xmin                            %
% beta_2 u + lambda_2 du/dx = g_2 for x = xmax                            %
%-------------------------------------------------------------------------%
% Discrete system :                                                       %
% u(x) = sum_i alpha_i R_ij                                               %
% (DE) lhs alpha = rhs                                                    %
%-------------------------------------------------------------------------%
% This code is distributed under the GNU gpl-3.0 license.                 %
%-------------------------------------------------------------------------%
% Input Parameters :                                                      %
% => Anonymous Function, a(x), b(x) : system coefficients                 %
% => Anonymous Function, f(x) : source term                               %
% => Integer, nx : number of space points                                 %
% => Real(nx), x : space points                                           %
% => Real, c : radial basis parameter                                     %
% Output Parameters :                                                     %
% => Real(nx), u : vector solution                                        %
% => Real(nx,nx), lhs : matrix of the discrete system                     %
% => Real(nx,1), rhs : vector of the source term and BC                   %
%-------------------------------------------------------------------------%
% Radial basis function
R = @(x,xi,c) ((x-xi).^2 + c.^2).^0.5;
% Solution nodes p = [x-coordinate, boundary nodes (=1)]
p = [x(:) (x(:)==x(1))+2*(x(:)==x(end))]; % In vector form
% In loop form
% for i=1:nx+1
%     p(i,1) = x(i); % x-coordinate
%     % Define boudary nodes
%     if (p(i,1) == x(1))
%         p(i,2)=1;
%     elseif(p(i,1) == x(end))
%         p(i,2)=2;
%     else
%         p(i,2)=0;
%     end
% end
% Construction of lhs, rhs and RB matrix
X = x'*ones(size(x));
Rij = R(X,X',c);
dRx = (X - X')./R(X,X',c); % dR/dx
d2Rx = c^2./(R(X,X',c)).^3; % d^2R/dx^2
lhs = diag(p(:,2)==1)*(beta(1)*Rij + lambda(1)*(X - X')./R(X,X',c)) + ...
      diag(p(:,2)==2)*(beta(2)*Rij + lambda(2)*(X - X')./R(X,X',c)) + ...
      diag(p(:,2)==0)*(a1(X).*d2Rx + a2(X).*dRx + a3(X).*Rij);
rhs = (p(:,2)==1)*(g(1)) + ...
      (p(:,2)==2)*(g(2)) + ...
      (p(:,2)==0).*(f(x'));
% In loop form
% lhs = zeros(nx,nx); rhs = zeros(nx,1); % lhs and rhs initialization
% for i=1:nx+1
%     for j=1:nx+1
%         Rij(i,j) = R(x(i),x(j),c); % Radial basis matrix
%         if(p(i,2) == 1) % Boundary nodes
%             lhs(i,j) = beta(1)*Rij(i,j) + ...
%                  lambda(1)*(x(i) - x(j))./R(x(i),x(j),c); % BC
%             rhs(i,1) = g(1);
%         elseif(p(i,2) == 2) % Boundary nodes
%             lhs(i,j) = beta(2)*Rij(i,j) + ...
%                  lambda(2)*(x(i) - x(j))./R(x(i),x(j),c); % BC
%             rhs(i,1) = g(2);
%         else % Inside nodes
%             dRx = (x(i) - x(j))./R(x(i),x(j),c); % dR/dx
%             d2Rx = c^2/(R(x(i),x(j),c)).^3; % d^2R/dx^2
%             lhs(i,j) = a1(x(i))*d2Rx + a2(x(i))*dRx + a3(x(i))*Rij(i,j);
%             rhs(i,1) = f(x(i));
%         end
%     end
% end
% Solution of discrete problem
alpha = lhs^-1*rhs;
% Solution construction
u = Rij*alpha;
