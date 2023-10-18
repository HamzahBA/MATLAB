clear all
close all
%------------------------------------------------------------------------
%% Physical, Mesh & Structure data
load PMData.mat
%------------------------------------------------------------------------
%% Artery geometry
load('IniStr.mat','Wup','Wlow')
% upper side
Reup = Wup(xup);
Rup = R0+[0*find(x<dup) Reup 0*find(x>dup+L0up)]; %geometry values
% lower side
Relow = Wlow(xlow);
Rlow = R0+[0*find(x<dlow) Relow 0*find(x>dlow+L0low)]; %geometry values
% sides splitting
yup = [zeros(size(find(y<0))) y(y>=0)];
ylow = [y(y<0) zeros(size(find(y>=0)))];
X = x'*ones(size(y)); Y = (Rup')*yup + (Rlow')*ylow; % meshgrids
yup = yup(2:end-1); ylow = ylow(2:end-1);
nsteps = 10;  % number of steps with graphic output
%------------------------------------------------------------------------
%% Initial & boundary conditions
U = zeros(Nx-2,Ny-2); V = zeros(Nx-2,Ny-2);
% boundary conditions
uwall = 0*x'; vwall = 0*x(2:end-1)';
u0 = 1; uin = u0*(1-y.^2); vin = y*0;
%------------------------------------------------------------------------
%% NS Operators
avg = @(M) (M(2:end,:)+M(1:end-1,:))/2; % values at middle points
%------------------------------------------------------------------------
disp('initialization')
%------------------------------------------------------------------------
% boudary coefficients 
uupx = 0; udownx = 1;
vupx = 0; vdownx = 1;
pupx = 1; pdownx = 0;
qupx = 0; qdownx = 0;
uupy = 0; udowny = 0;
vupy = 0; vdowny = 0;
pupy = 1; pdowny = 1;
qupy = 0; qdowny = 0;
% central approximation of Laplacian operator
D2 = @(n,h,up,down) ...
    spdiags([[1 -2+up 0];ones(n-2,1)*[1 -2 1];[0 -2+down 1]],-1:1,n,n)/h^2;
% central approximation of first derivative
D1 = @(n,h,up,down) ...
    spdiags([[-1 -up 1];ones(n-2,1)*[-1 0 1];[0 down 1]],-1:1,n,n)/2/h;
% central approximation of normalized Laplacian operator
% d^2/dX^2
LapUbc = [uin(2:end-1);zeros(Nx-3,Ny-2)]/dx^2;
uLapx1 = kron(speye(Ny-2),D2(Nx-2,dx,uupx,udownx));
vLapx1 = kron(speye(Ny-2),D2(Nx-2,dx,vupx,vdownx));
pLapx1 = kron(speye(Ny-2),D2(Nx-2,dx,pupx,pdownx));
qLapx1 = kron(speye(Ny-2),D2(Nx-2,dx,qupx,qdownx));
% -Y*dR/dX*d(1/R)/dX*d/dY
Lapx21 = @(n,y,up,down) spdiags(-y',0,n,n)*D1(Ny-2,dy,up,down);
Lapx22 = @(n,R) spdiags(diff(avg(R')).*diff(avg(1./(R')))/dx^2,0,n,n);
uLapx2 = kron(Lapx21(Ny-2,yup,uupy,udowny),Lapx22(Nx-2,Rup)) + ...
    kron(Lapx21(Ny-2,ylow,uupy,udowny),Lapx22(Nx-2,Rlow));
vLapx2 = kron(Lapx21(Ny-2,yup,vupy,vdowny),Lapx22(Nx-2,Rup)) + ...
    kron(Lapx21(Ny-2,ylow,vupy,vdowny),Lapx22(Nx-2,Rlow));
pLapx2 = kron(Lapx21(Ny-2,yup,pupy,pdowny),Lapx22(Nx-2,Rup)) + ...
    kron(Lapx21(Ny-2,ylow,pupy,pdowny),Lapx22(Nx-2,Rlow));
qLapx2 = kron(Lapx21(Ny-2,yup,qupy,qdowny),Lapx22(Nx-2,Rup)) + ...
    kron(Lapx21(Ny-2,ylow,qupy,qdowny),Lapx22(Nx-2,Rlow));
% -Y/R*d^2R/dX^2*d/dY
Lapx31 = @(n,y,up,down) spdiags(-y',0,n,n)*D1(Ny-2,dy,up,down);
Lapx32 = @(n,R) spdiags(diff(R,2)'./(R(2:end-1)')/dx^2,0,n,n);
uLapx3 = kron(Lapx31(Ny-2,yup,uupy,udowny),Lapx32(Nx-2,Rup)) + ...
    kron(Lapx31(Ny-2,ylow,uupy,udowny),Lapx32(Nx-2,Rlow));
vLapx3 = kron(Lapx31(Ny-2,yup,vupy,vdowny),Lapx32(Nx-2,Rup)) + ...
    kron(Lapx31(Ny-2,ylow,vupy,vdowny),Lapx32(Nx-2,Rlow));
pLapx3 = kron(Lapx31(Ny-2,yup,pupy,pdowny),Lapx32(Nx-2,Rup)) + ...
    kron(Lapx31(Ny-2,ylow,pupy,pdowny),Lapx32(Nx-2,Rlow));
qLapx3 = kron(Lapx31(Ny-2,yup,qupy,qdowny),Lapx32(Nx-2,Rup)) + ...
    kron(Lapx31(Ny-2,ylow,qupy,qdowny),Lapx32(Nx-2,Rlow));
% -2*Y/R*dR/dX*d^2/dX/dY
LapUbc = LapUbc + [(-yup*(Rup(3)-Rup(1))/2/dx/(Rup(2)) ...
    -ylow*(Rlow(3)-Rlow(1))/2/dx/(Rlow(2))).* ...
    (-diff(avg(uin')')/dy);zeros(Nx-3,Ny-2)]/dx;
Lapx41 = @(n,y,up,down) spdiags(-y',0,n,n)*D1(Ny-2,dy,up,down);
Lapx42 = @(n,R,up,down) ...
    spdiags(2*diff(avg(R'))./(R(2:end-1)')/dx,0,n,n)*D1(Nx-2,dx,up,down);
uLapx4 = kron(Lapx41(Ny-2,yup,uupy,udowny),Lapx42(Nx-2,Rup,uupx,udownx)) + ...
    kron(Lapx41(Ny-2,ylow,uupy,udowny),Lapx42(Nx-2,Rlow,uupx,udownx));
vLapx4 = kron(Lapx41(Ny-2,yup,vupy,vdowny),Lapx42(Nx-2,Rup,vupx,vdownx)) + ...
    kron(Lapx41(Ny-2,ylow,vupy,vdowny),Lapx42(Nx-2,Rlow,vupx,vdownx));
pLapx4 = kron(Lapx41(Ny-2,yup,pupy,pdowny),Lapx42(Nx-2,Rup,pupx,pdownx)) + ...
    kron(Lapx41(Ny-2,ylow,pupy,pdowny),Lapx42(Nx-2,Rlow,pupx,pdownx));
qLapx4 = kron(Lapx41(Ny-2,yup,qupy,qdowny),Lapx42(Nx-2,Rup,qupx,qdownx)) + ...
    kron(Lapx41(Ny-2,ylow,qupy,qdowny),Lapx42(Nx-2,Rlow,qupx,qdownx));
% Y/R^2*(dR/dX)^2*d/dY
Lapx51 = @(n,y,up,down) spdiags(y',0,n,n)*D1(Ny-2,dy,up,down);
Lapx52 = @(n,R,up,down) ...
    spdiags(diff(avg(R')).^2./(R(2:end-1)'.^2)/dx^2,0,n,n);
uLapx5 = kron(Lapx51(Ny-2,yup,uupy,udowny),Lapx52(Nx-2,Rup)) + ...
    kron(Lapx51(Ny-2,ylow,uupy,udowny),Lapx52(Nx-2,Rlow));
vLapx5 = kron(Lapx51(Ny-2,yup,vupy,vdowny),Lapx52(Nx-2,Rup)) + ...
    kron(Lapx51(Ny-2,ylow,vupy,vdowny),Lapx52(Nx-2,Rlow));
pLapx5 = kron(Lapx51(Ny-2,yup,pupy,pdowny),Lapx52(Nx-2,Rup)) + ...
    kron(Lapx51(Ny-2,ylow,pupy,pdowny),Lapx52(Nx-2,Rlow));
qLapx5 = kron(Lapx51(Ny-2,yup,qupy,qdowny),Lapx52(Nx-2,Rup)) + ...
    kron(Lapx51(Ny-2,ylow,qupy,qdowny),Lapx52(Nx-2,Rlow));
% Y^2/R^2*(dR/dX)^2*d^2/dY^2
Lapx61 = @(n,y,up,down) spdiags(y'.^2,0,n,n)*D2(Ny-2,dy,up,down);
Lapx62 = @(n,R) spdiags(diff(avg(R')).^2./(R(2:end-1)'.^2)/dx^2,0,n,n);
uLapx6 = kron(Lapx61(Ny-2,yup,uupy,udowny),Lapx62(Nx-2,Rup)) + ...
    kron(Lapx61(Ny-2,ylow,uupy,udowny),Lapx62(Nx-2,Rlow));
vLapx6 = kron(Lapx61(Ny-2,yup,vupy,vdowny),Lapx62(Nx-2,Rup)) + ...
    kron(Lapx61(Ny-2,ylow,vupy,vdowny),Lapx62(Nx-2,Rlow));
pLapx6 = kron(Lapx61(Ny-2,yup,pupy,pdowny),Lapx62(Nx-2,Rup)) + ...
    kron(Lapx61(Ny-2,ylow,pupy,pdowny),Lapx62(Nx-2,Rlow));
qLapx6 = kron(Lapx61(Ny-2,yup,qupy,qdowny),Lapx62(Nx-2,Rup)) + ...
    kron(Lapx61(Ny-2,ylow,qupy,qdowny),Lapx62(Nx-2,Rlow));
% d^2/dx^2 = d^2/dX^2 + Y^2/R^2*(dR/dX)^2*d^2/dY^2
%            - Y/R*d^2R/dX^2*d/dY - Y/R*dR/dX*d^2/dX/dY
%            + Y/R^2*(dR/dX)^2*d/dY + Y^2/R^2*(dR/dX)^2*d^2/dY^2
uLapx = uLapx1 + uLapx2 + uLapx3 + uLapx4 + uLapx5 + uLapx6;
vLapx = vLapx1 + vLapx2 + vLapx3 + vLapx4 + vLapx5 + vLapx6;
pLapx = pLapx1 + pLapx2 + pLapx3 + pLapx4 + pLapx5 + pLapx6;
qLapx = qLapx1 + qLapx2 + qLapx3 + qLapx4 + qLapx5 + qLapx6;
% d^2/dy^2 = 1/R^2*d^2/dY^2
Lapy = @(n,R) spdiags(1./(R(2:end-1)'.^2),0,n,n);
uLapy = kron(spdiags((y(2:end-1)'>=0),0,Ny-2,Ny-2)*D2(Ny-2,dy,uupy,udowny),Lapy(Nx-2,Rup)) + ...
    kron(spdiags((y(2:end-1)'<0),0,Ny-2,Ny-2)*D2(Ny-2,dy,uupy,udowny),Lapy(Nx-2,Rlow));
ULap = uLapx + uLapy;
    
VLap = vLapx + ...
    kron(spdiags((y(2:end-1)'>=0),0,Ny-2,Ny-2)*D2(Ny-2,dy,vupy,vdowny),Lapy(Nx-2,Rup)) + ...
    kron(spdiags((y(2:end-1)'<0),0,Ny-2,Ny-2)*D2(Ny-2,dy,vupy,vdowny),Lapy(Nx-2,Rlow));
PLap = pLapx + ...
    kron(spdiags((y(2:end-1)'>=0),0,Ny-2,Ny-2)*D2(Ny-2,dy,pupy,pdowny),Lapy(Nx-2,Rup)) + ...
    kron(spdiags((y(2:end-1)'<0),0,Ny-2,Ny-2)*D2(Ny-2,dy,pupy,pdowny),Lapy(Nx-2,Rlow));
QLap = qLapx + ...
    kron(spdiags((y(2:end-1)'>=0),0,Ny-2,Ny-2)*D2(Ny-2,dy,qupy,qdowny),Lapy(Nx-2,Rup)) + ...
    kron(spdiags((y(2:end-1)'<0),0,Ny-2,Ny-2)*D2(Ny-2,dy,qupy,qdowny),Lapy(Nx-2,Rlow));
% d/dX
DxUbc = [-uin(2:end-1);zeros(Nx-3,Ny-2)]/2/dx;
uDx1 = kron(speye(Ny-2),D1(Nx-2,dx,uupx,udownx));
vDx1 = kron(speye(Ny-2),D1(Nx-2,dx,vupx,vdownx));
% -Y/R*dR/dX*d/dY
Dx21 = @(n,y,up,down) spdiags(-y',0,n,n)*D1(Ny-2,dy,up,down);
Dx22 = @(n,R) ...
    spdiags(diff(avg(R'))./(R(2:end-1)')/dx,0,n,n);
uDx2 = kron(Dx21(Ny-2,yup,uupy,udowny),Dx22(Nx-2,Rup)) + ...
    kron(Dx21(Ny-2,ylow,uupy,udowny),Dx22(Nx-2,Rlow));
vDx2 = kron(Dx21(Ny-2,yup,vupy,vdowny),Dx22(Nx-2,Rup)) + ...
    kron(Dx21(Ny-2,ylow,vupy,vdowny),Dx22(Nx-2,Rlow));
% d/dx = d/dX - Y/R*dR/dX*d/dY
% d/dy = 1/R*d/dY
uDx = uDx1 + uDx2;
vDx = vDx1 + vDx2;
% d/dy = 1/R*d/dY
Dy = @(n,R) spdiags(1./(R(2:end-1)'),0,n,n);
uDy = kron(spdiags((y(2:end-1)'>=0),0,Ny-2,Ny-2)*D1(Ny-2,dy,uupy,udowny),Dy(Nx-2,Rup)) + ...
    kron(spdiags((y(2:end-1)'<0),0,Ny-2,Ny-2)*D1(Ny-2,dy,uupy,udowny),Dy(Nx-2,Rlow));
vDy = kron(spdiags((y(2:end-1)'>=0),0,Ny-2,Ny-2)*D1(Ny-2,dy,vupy,vdowny),Dy(Nx-2,Rup)) + ...
    kron(spdiags((y(2:end-1)'<0),0,Ny-2,Ny-2)*D1(Ny-2,dy,vupy,vdowny),Dy(Nx-2,Rlow));
% normalized derivatives
diffx = @(u,v) diff(avg(u))/dx - (((1./(Rup(2:end-1)')).* ...
    diff(avg(Rup'))/dx)*yup).*diff(avg(v'))'/dy - (((1./(Rlow(2:end-1)')).* ...
    diff(avg(Rlow'))/dx)*ylow).*diff(avg(v'))'/dy;
diffy = @(u) (1./Rup(2:end-1)'*((y(2:end-1)>=0))).*(diff(avg(u'))'/dy) + ...
    (1./Rlow(2:end-1)'*((y(2:end-1)<0))).*(diff(avg(u'))'/dy);
%------------------------------------------------------------------------
%% Time iterations
Psol=[];Qsol=[];
for i=1:Nt
    % treat nonlinear terms
    % semi-implicit convection/ implicit diffusion
    UConv = spdiags(U(:),0,(Nx-2)*(Ny-2),(Nx-2)*(Ny-2))*uDx ...
            + spdiags(V(:),0,(Nx-2)*(Ny-2),(Nx-2)*(Ny-2))*uDy;
    ConvUbc = spdiags(U(:),0,(Nx-2)*(Ny-2),(Nx-2)*(Ny-2))*DxUbc(:);
    Urhs = (speye((Nx-2)*(Ny-2)) - dt/Re*ULap + dt*UConv) ...
        \(U(:)+ dt/Re*LapUbc(:) - dt*ConvUbc);
    VConv = spdiags(U(:),0,(Nx-2)*(Ny-2),(Nx-2)*(Ny-2))*vDx ...
            + spdiags(V(:),0,(Nx-2)*(Ny-2),(Nx-2)*(Ny-2))*vDy;
    Vrhs = (speye((Nx-2)*(Ny-2)) - dt/Re*VLap + dt*VConv)\(V(:));
    U = reshape(Urhs,Nx-2,Ny-2);
    V = reshape(Vrhs,Nx-2,Ny-2);
    % pressure correction
    Ue = [uin(2:end-1);U;U(end,:)]; Ue = [uwall Ue uwall];
    Ve = [vwall V vwall]; Ve = [vin;Ve;Ve(end,:)];
    Ux = diffx(Ue(:,2:end-1),Ue(2:end-1,:)); Uy = diffy(Ue(2:end-1,:));
    Vx = diffx(Ve(:,2:end-1),Ve(2:end-1,:)); Vy = diffy(Ve(2:end-1,:));
    rhs = PLap\reshape((Ux+Vy)/dt,[],1);
    P = reshape(rhs,Nx-2,Ny-2);
    Pe = [P(1,:);P;0*y(2:end-1)]; Pe = [Pe(:,1) Pe Pe(:,end)];
    Px = diffx(Pe(:,2:end-1),Pe(2:end-1,:)); Py = diffy(Pe(2:end-1,:));
    U = U-dt*Px;
    V = V-dt*Py;
    % adding boundary conditions
    Papp = [P(1,:);P;0*y(2:end-1)]; Papp = [Papp(:,1) Papp Papp(:,end)];
    Uapp = [uin(2:end-1);U;U(end,:)]; Uapp = [uwall Uapp uwall];
    Vapp = [vwall V vwall]; Vapp = [vin;Vapp;Vapp(end,:)];
    % stream function
    qin = -u0*R0*(y-y.^3/3) - 2/3*u0*R0;
    qup = -4/3*u0*R0 + 0*x';%trapz(R0*y,-Uapp(end,:)) - 2/3*u0*R0 + 0*x';%[-4/3*u0*R0  + 0*x'];
    for j=2:Ny
        qout(j) = trapz(R0*y(1:j),-R0*Uapp(end,1:j)) - 2/3*u0*R0;
    end
    Qbc = [qin(2:end-1);zeros(Nx-4,Ny-2);qout(2:end-1)]/dx^2;
    Qbc = Qbc + [zeros(Nx-2,Ny-3) ...
        (-y(end-1).*diff(avg(1./Rup'))/dx.*diff(avg(Rup'))/dx.*qup(2:end-1))]/2/dy;
    Qbc = Qbc + [zeros(Nx-2,Ny-3) ...
        (-y(end-1)./Rup(2:end-1)'.*diff(Rup,2)'/dx^2.*qup(2:end-1))]/2/dy;
    Qbc = Qbc + [zeros(Nx-2,Ny-3) ...
        (-2*y(end-1)./Rup(2:end-1)'.*diff(avg(Rup'))/dx.*diff(avg(qup))/dx)]/2/dy;
    Qbc = Qbc + [(-2*yup./Rup(2).*((Rup(3)-Rup(1))/2/dx) ...
        -2*ylow./Rlow(2).*((Rlow(3)-Rlow(1))/2/dx)).*(-diff(avg(qin')')/dy)...
        ;zeros(Nx-4,Ny-2); ...
        (-2*yup./Rup(Nx-1).*((Rup(Nx)-Rup(Nx-2))/2/dx) ...
        -2*ylow./Rlow(Nx-1).*((Rlow(Nx)-Rlow(Nx-2))/2/dx)).*(diff(avg(qout')')/dy)]/2/dx;
    Qbc = Qbc + [zeros(Nx-2,Ny-3) ...
        (y(end-1)./Rup(2:end-1)'.^2.*(diff(avg(Rup'))/dx).^2.*qup(2:end-1))]/2/dy;
    Qbc = Qbc + [zeros(Nx-2,Ny-3) ...
        (y(end-1).^2./Rup(2:end-1)'.^2.*(diff(avg(Rup'))/dx).^2.*qup(2:end-1))]/dy^2;
    Qbc = Qbc + [zeros(Nx-2,Ny-3) (1./Rup(2:end-1)'.^2.*qup(2:end-1))]/dy^2;
    Uy = diffy(Uapp(2:end-1,:));
    Vx = diffx(Vapp(:,2:end-1),Vapp(2:end-1,:));
    rhs = -QLap\reshape(Uy-Vx+Qbc,[],1);
    Q = reshape(rhs,Nx-2,Ny-2);
    Q = [qin(2:end-1);Q;qout(2:end-1)];
    Q = [0*x' Q qup]; Qapp = Q;
    % saving values from selected iterations
    if mod((i-1),Nt/nsteps) == 1, fprintf('.'), end
%     if mod((i-1),nsteps) == 1
%         Psol = [Psol Papp(:)];
%         Qsol = [Qsol Q(:)];
%     end
end
%------------------------------------------------------------------------
%% Saving results
save('InitPress.mat','dt','X','Y','Re','Uapp','Vapp','Papp','Qapp')