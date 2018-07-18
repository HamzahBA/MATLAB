function [ u ] = RD_MeshlessRBF1d (  c, x, nx )

%*****************************************************************************80
%
%% RBF meshless code for the 1D reaction-diffusion equation
%
%  Discussion:
%
%    The solution is defined over the unit square [0,1]x[0,1].
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    23 January 2015
%
%  Author:
%
%    John Burkardt
%
%  Reference:
%
%    Junping Wang, Yanqiu Wang, Xiu Ye,
%    A robust numerical method for Stokes equations based on divergence-free
%    H(div) finite element methods,
%    SIAM Journal on Scientific Computing,
%    Volume 31, Number 4, 2009, pages 2784-2802.
%
%  Parameters:
%
%    Input, integer N, the number of evaluation points.
%
%    Input, real X(N), Y(N), the coordinates of the points.
%
%    Output, real U(N), V(N), P(N), the velocity components and
%    pressure at each of the points.

clear all
close all

xmin = 0; xmax = 1;
nx = 100; dx = (xmax-xmin)/nx;

x = xmin:dx:xmax; X = x';
p = [X(:) (X(:)==xmin | X(:)==xmax)];

c = 0.07;

f = @(x,y) pi^2*sin(pi*x);

R = @(x,xi,c) ((x-xi).^2 + c.^2).^0.5;

for i=1:nx+1
    for j=1:nx+1
        Rij(i,j) = R(X(i),X(j),c);
        switch p(i,2)
        case 1
            A(i,j) = Rij(i,j);
            F(i,1) = 0;
        case 0
            d2Rx = c^2/(R(X(i),X(j),c)).^3;
            A(i,j) = -d2Rx;
            F(i,1) = f(X(i));
        end
    end
end

alpha = A^-1*F;

u = Rij*alpha;

uex = sin(pi*X);

plot(x,u,x,uex,'*r')

