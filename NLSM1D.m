clear all
close all
%------------------------------------------------------------------------
%% Physical, Mesh & Structure data
load PMData.mat
%------------------------------------------------------------------------
%% Steady Structure

Wini = @(P) [(lambda_s/2 + mu_s)*P^2;0;P];
F = @(x,W,f) [f(x)/thick*W(3)/sqrt(1+W(3)^2); ...
    W(3); ...
    1/(W(1) + mu_s)*(-f(x)/thick/sqrt(1+W(3)^2))];
opt=odeset('RelTol',1e-5);

% Upper side
fup = @(x) P0int - Pext;
p1up = 0; p2up = p1up + 0.1;
W0p1up = Wini(p1up); W0p2up = Wini(p2up);

[xxup,Wp1up] = ode45(@(x,W) F(x,W,fup), xxup, W0p1up);
[xxup,Wp2up] = ode45(@(x,W) F(x,W,fup), xxup, W0p2up);

eps = 10^-15;n = 10000;
err = 1;k = 0;
while(err > eps & k<n)
    
    fp1up = Wp1up(end,2); fp2up = Wp2up(end,2);
    p3up = p2up - (p2up-p1up)/(fp2up-fp1up)*fp2up;
    
    W0p3up = Wini(p3up);
    [xxup,Wp3up] = ode45(@(x,W) F(x,W,fup), xxup, W0p3up);
    Wp3up;
    
    err = norm(p3up-p2up);
    p1up = p2up; p2up = p3up;
    Wp1up = Wp2up; Wp2up = Wp3up;
    k = k + 1;
end

Wup = @(val) interp1(xxup',Wp3up(:,2),val);

% Lower side
flow = @(x) P0int - Pext;
p1low = 0; p2low = p1low + 0.1;
W0p1low = Wini(p1low); W0p2low = Wini(p2low);

[xxlow,Wp1low] = ode45(@(x,W) F(x,W,flow), xxlow, W0p1low);
[xxlow,Wp2low] = ode45(@(x,W) F(x,W,flow), xxlow, W0p2low);

err = 1;k = 0;
while(err > eps & k<n)
    
    fp1low = Wp1low(end,2); fp2low = Wp2low(end,2);
    p3low = p2low - (p2low-p1low)/(fp2low-fp1low)*fp2low;
    
    W0p3low = Wini(p3low);
    [xxlow,Wp3low] = ode45(@(x,W) F(x,W,flow), xxlow, W0p3low);
    
    err = norm(p3low-p2low);
    p1low = p2low; p2low = p3low;
    Wp1low = Wp2low; Wp2low = Wp3low;
    k = k + 1;
end

Wlow = @(val) interp1(xxlow',Wp3low(:,2),val);
%------------------------------------------------------------------------
% figure
% pii = (lambda_s/2 + mu_s)*Wp3low(:,3).^2;
% plot(xlow,1./(pii+mu_s)*flow(x)/thick)
% figure
% plot(xlow,Wp3low(:,3))
figure
plot(xup,Wup(xup))
figure
plot(xlow,Wlow(xlow))
%------------------------------------------------------------------------
%% Saving results
save('IniStr.mat','Wup','Wlow')