clear all
close all
%% Physical data
mu = 3.5*10^-3; rho = 1.06*10^3; D0 = 8*10^-2*10; Re = 200; meanU = mu*Re/rho/D0;
%Structural data
rho_s = 1*10^-3/10^-6; thick = 0.002; mu_s = 1.1/10^-6; lambda_s = 1.7/10^-6;
P0int = 13000; Pext = 5000;
% Non-dimensional parameters
% P0int = P0int/rho/meanU^2/10000; Pext = Pext/rho/meanU^2/10000;
% lambda_s = lambda_s/rho/meanU^2/100000; mu_s = mu_s/rho/meanU^2/100000;
% thick = thick/D0*10; rho_s = rho_s/rho;

P0int = 15; Pext = 5; rho_s = 1;
lambda_s = 1000; mu_s = 1000;
thick = 0.05; rho_s = 1;
%% Mesh configurations
%------------------------------------------------------------------------
t0 = 0; T = 50;    % final time
R0 = 0.5; L = 10;
Nx = 101; Ny = 51;      % number of x and y gridpoints
Nt = 400;
x = linspace(0,L,Nx); dx = L/(Nx-1);
y = linspace(-1,1,Ny); dy = 2/(Ny-1);
t = linspace(t0,T,Nt+1); dt = (T-t0)/Nt;
%------------------------------------------------------------------------
%% Initial artery geometry
% upper side
dup = 0.1*L; L0up = 0.4*L; L0 = L0up; % first deformable part
xeup = x(find(x>=dup & x<=dup+L0up));
xup = xeup - dup;
xxup = 0:0.1:L0up;
% lower side
dlow = 0.25*L; L0low = 0.3*L; % second deformable part
xelow = x(find(x>=dlow & x<=dlow+L0low));
xlow = xelow - dlow;
xxlow = 0:0.1:L0low;
%------------------------------------------------------------------------
save PMData