clear all
close all
%------------------------------------------------------------------------
%% Physical, Mesh & Structure data
load PMData.mat
%------------------------------------------------------------------------
%% Steady Structure
xx = xup;

Wini = @(P) [(lambda_s/2 + mu_s)*P^2;0;P];
F = @(x,W,f) [f(x)/thick*W(3)/sqrt(1+W(3)^2); ...
    W(3); ...
    1/(W(1) + mu_s)*(-f(x)/thick*sqrt(1+W(3)^2))];
opt=odeset('RelTol',1e-5);

% Upper side
f = @(x) P0int - Pext;
p1 = 0; p2 = p1 + 0.1;
W0p1 = Wini(p1); W0p2 = Wini(p2);

[xx,Wp1] = RK4(@(x,W) F(x,W,f), xx, W0p1);
[xx,Wp2] = RK4(@(x,W) F(x,W,f), xx, W0p2);
% [xx,Wp1] = ode45(@(x,W) F(x,W,f), xx, W0p1,opt);
% [xx,Wp2] = ode45(@(x,W) F(x,W,f), xx, W0p2,opt);

eps = 10^-15;n = 10000;
err = 1;k = 0;

while(err > eps & k<n)
    
    fp1 = Wp1(end,2); fp2 = Wp2(end,2);
    p3 = p2 - (p2-p1)/(fp2-fp1)*fp2;
    
    W0p3 = Wini(p3);
    [xx,Wp3] = RK4(@(x,W) F(x,W,f), xx, W0p3);
%     [xx,Wp3] = ode45(@(x,W) F(x,W,f), xx, W0p3,opt);
    Wp3;
    
    err = norm(p3-p2);
    p1 = p2; p2 = p3;
    Wp1 = Wp2; Wp2 = Wp3;
    k = k + 1;
end

W = @(val) interp1(xx',Wp3(:,2),val);

plot(xx,W(xx))