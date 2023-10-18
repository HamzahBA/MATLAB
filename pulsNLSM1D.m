clear all
close all
%------------------------------------------------------------------------
%% Physical, Mesh & Structure data
load PMData.mat
load('IniStr.mat','Wup','Wlow')
Wupold = Wup; Wlowold = Wlow;
%------------------------------------------------------------------------
%% Steady Structure
P = @(t) 20*sin(4*pi*t/50);
for i=1:Nt
    Pup = @(val) interp1(x',P0int + P(t(i+1)) + 0*x',val);
    Plow = @(val) interp1(x',P0int + P(t(i+1)) + 0*x',val);
    
    fup = @(x) Pup(x+dup) - Pext;
    flow = @(x) Plow(x+dlow) - Pext;
    
    p1up = 0; p2up = p1up + 0.1;
    p1low = 0; p2low = p1low + 0.1;
    
    Wini = @(P) [(lambda_s/2 + mu_s)*P^2;0;P];
    W0p1up = Wini(p1up); W0p2up = Wini(p2up);
    W0p1low = Wini(p1low); W0p2low = Wini(p2low);
    
    F = @(x,W,Wcurr,Wold,f) [f(x)/thick*W(3)/sqrt(1+W(3)^2); ...
        W(3); ...
        1/(W(1) + mu_s)*(rho_s*(W(2)-2*Wcurr(x)+Wold(x))/dt^2-f(x)/thick*sqrt(1+W(3)^2))];
    
    opt=odeset('RelTol',1e-5);
    
    [xxup,Wp1up] = ode45(@(x,W) F(x,W,Wup,Wupold,fup), xup, W0p1up,opt);
    [xxup,Wp2up] = ode45(@(x,W) F(x,W,Wup,Wupold,fup), xup, W0p2up,opt);
    
    [xxlow,Wp1low] = ode45(@(x,W) F(x,W,Wlow,Wlowold,flow), xlow, W0p1low,opt);
    [xxlow,Wp2low] = ode45(@(x,W) F(x,W,Wlow,Wlowold,flow), xlow, W0p2low,opt);
    
    eps = 10^-15;n = 10000;
    err = 1;k = 0;
    while(err > eps & k<n)
        
        fp1up = Wp1up(end,2); fp2up = Wp2up(end,2);
        p3up = p2up - (p2up-p1up)/(fp2up-fp1up)*fp2up;
        
        W0p3up = Wini(p3up);
        [xxup,Wp3up] = ode45(@(x,W) F(x,W,Wup,Wupold,fup), xup, W0p3up,opt);
        
        err = norm(p3up-p2up);
        p1up = p2up; p2up = p3up;
        Wp1up = Wp2up; Wp2up = Wp3up;
        
        k = k + 1;
    end
    
    SM = F(xup(end),Wp3up(end,:),Wup,Wupold,fup); SM(3)
    
    err = 1;k = 0;
    while(err > eps & k<n)
        
        fp1low = Wp1low(end,2); fp2low = Wp2low(end,2);
        p3low = p2low - (p2low-p1low)/(fp2low-fp1low)*fp2low;
        
        W0p3low = Wini(p3low);
        [xxlow,Wp3low] = ode45(@(x,W) F(x,W,Wlow,Wlowold,flow), xlow, W0p3low,opt);
        
        err = norm(p3low-p2low);
        p1low = p2low; p2low = p3low;
        Wp1low = Wp2low; Wp2low = Wp3low;
        
        k = k + 1;
    end
    
    % upper side
    Wupold = Wup;
    Wup = @(val) interp1(xup',Wp3up(:,2),val);
    % lower side
    Wlowold = Wlow;
    Wlow = @(val) interp1(xlow',Wp3low(:,2),val);
    
%     figure
    plot(xup,Wup(xup))
%     figure
%     plot(xlow,Wlow(xlow))
    drawnow
end
%------------------------------------------------------------------------
% figure
% pii = (lambda_s/2 + mu_s)*Wp3low(:,3).^2;
% plot(xlow,1./(pii+mu_s)*flow(x)/thick)
% figure
% plot(xlow,Wp3low(:,3))
% figure
% plot(xup,Wup(xup))
% figure
% plot(xlow,Wlow(xlow))
%------------------------------------------------------------------------
%% Saving results
save('IniStr.mat','Wup','Wlow')