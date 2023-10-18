clear all
close all

load('SymTwoAneurysmDef.mat','u0pulse','dt','Nx','Ny','Re','Psol','Qsol','Xsol','Ysol')
% load('NSymTwoAneurysmDef.mat','dt','X','Y','Re','Psol','Qsol')
% load('NSymTwoAneurysmDef1.mat','dt','X','Y','dlow','Psol','Qsol')
% X1 = X; Y1 = Y; dlow1 = dlow; Psol1 = Psol; Qsol1 = Qsol;
% load('NSymTwoAneurysmDef2.mat','dt','X','Y','dlow','Psol','Qsol')
% X2 = X; Y2 = Y; dlow2 = dlow; Psol2 = Psol; Qsol2 = Qsol;
% load('NSymTwoAneurysmDef3.mat','dt','X','Y','dlow','Psol','Qsol')
% X3 = X; Y3 = Y; dlow3 = dlow; Psol3 = Psol; Qsol3 = Qsol;
% [Nx,Ny] = size(X);

u0 = u0pulse;
Nt = length(Psol(1,:));
xmin = 0.5; xmax = 6;
ymin = min(Ysol(:)); ymax = max(Ysol(:));
L = max(Xsol(:)); dx = (L)/(Nx-1);
Nx0 = round(xmin*(Nx-1)/L); Nx1 = round(xmax*(Nx+1)/L);
% X1 = X1(Nx0:Nx1,:); Y1 = Y1(Nx0:Nx1,:);
% X2 = X2(Nx0:Nx1,:); Y2 = Y2(Nx0:Nx1,:);
% X3 = X3(Nx0:Nx1,:); Y3 = Y3(Nx0:Nx1,:);
% figure
% t = 0:dt:50;
% hu0 = plot(t,u0,'r', 'LineWidth', 2); hold on;
% plot(t(round([18:2:26 31]*Nt/50)),u0(round([18:2:26 31]*Nt/50)),'ok' ...
%    , 'MarkerFaceColor', 'k', 'MarkerSize', 10);
% text(t(round([18:2:26 31]*Nt/50)),u0(round([18:2:26 31]*Nt/50)), ...
%     ['(t_1)';'(t_2)';'(t_3)';'(t_4)';'(t_5)';'(t_6)'], 'HorizontalAlignment', 'right' ...
%     , 'VerticalAlignment', 'bottom', 'FontSize', 30);
% hold off; axis([15 31.5 min(u0)-0.1 max(u0)+0.2])
% % axis([10 50 min(u0)-0.01 max(u0)+0.01])
% set(gca, 'fontsize', 15);

X = reshape(Xsol(:,1),Nx,Ny); Y = reshape(Ysol(:,1),Nx,Ny);
X = X(Nx0:Nx1,:); Y = Y(Nx0:Nx1,:);
iupmin = 0; minup = norm(Y(:,end),1);
iupmax = 0; maxup = norm(Y(:,end),1);
ilowmin = 0; minlow = norm(Y(:,1),1);
ilowmax = 0; maxlow = norm(Y(:,1),1);
k = 0;
for i=1:Nt%round([18:2:26 31]*Nt/50)%round([18:2:26 31]*Nt/50)%round(10*Nt/50)%round([10:5:35]*Nt/50)%round([5 10 14.7 25]*Nt/50)%round(10*Nt/50)%round([1 3 4 6 8 10 15 20 25]*Nt/50)%round(10*Nt/20)%Nt%round([1 3 4.5 6 7 12]*Nt/20)%length(Psol(1,:))%round(9*Nt/20)
    % visualization
    %     stream function
        k = k + 1;
    %     subplot(4,1,k);
    X = reshape(Xsol(:,i),Nx,Ny); Y = reshape(Ysol(:,i),Nx,Ny);
    X = X(Nx0:Nx1,:); Y = Y(Nx0:Nx1,:);
    P = reshape(Psol(:,i),Nx,Ny); P = P(Nx0:Nx1,:);
    Q = reshape(Qsol(:,i),Nx,Ny); Q = Q(Nx0:Nx1,:);
%         Q1 = reshape(Qsol1(:,i),Nx,Ny); Q1 = Q1(Nx0:Nx1,:);
%         Q2 = reshape(Qsol2(:,i),Nx,Ny); Q2 = Q2(Nx0:Nx1,:);
%         Q3 = reshape(Qsol3(:,i),Nx,Ny); Q3 = Q3(Nx0:Nx1,:);
%     [p1,p2]=contourf(X,Y,P,30);set(p2,'edgecolor','none'); hold on
%     p = sort(P(:)); caxis(p([1 end])); colorbar()
% %     contour(X,Y,Q,60,'-k');hold on
% %     plot(X(:,1),Y(:,1),'b',X(:,end),Y(:,end),'b')
%     hold off, axis([xmin xmax ymin-0.1 ymax+0.1])
%     set(gca, 'fontsize', 15);
    
    if(norm(Y(:,1),1) < minlow)
        minlow = norm(Y(:,1),1);
        ilowmin = i;
    end
    if(norm(Y(:,1),1) > maxlow)
        maxlow = norm(Y(:,1),1);
        ilowmax = i;
    end
    if(norm(Y(:,end),1) < minup)
        minup = norm(Y(:,end),1);
        iupmin = i;
    end
    if(norm(Y(:,end),1) > maxup)
        maxup = norm(Y(:,end),1);
        iupmax = i;
    end
    
%     drawnow
%     
%     h=gcf;
%     set(h,'renderer','opengl');
%     set(h, 'PaperPosition', [0 0 10 6])
%     set(h,'PaperSize', [10 6]);
%     
% %     print(h,strcat('PulsNSymTwoElaAneurysmT',num2str(k),'SL'),'-dpdf')
%     print(h,strcat('PulsNSymTwoElaAneurysmT',num2str(k),'P'),'-dpdf')
%     
%     pause(1)

%     h=gcf;
% %     set(h,'renderer','opengl');
%     set(h, 'PaperPosition', [0 0 10 6])
%     set(h,'PaperSize', [10 6]);
%     
%     saveas(h,strcat('Animation/PulsNSymTwoElaAneurysmT',num2str(k)),'fig')
%     
%     pause(1)
end

figure
%     k = k + 1;
    X = reshape(Xsol(:,1),Nx,Ny); x = X(:,1);
%     Y = reshape(Ysol(:,1),Nx,Ny); y = Y(:,1);
    Y = reshape(Ysol(:,1),Nx,Ny); y = Y(:,end);
    plot(x,y, 'r-', 'LineWidth', 1); hold on
%     Y = reshape(Ysol(:,ilowmin),Nx,Ny); y = Y(:,1);
    Y = reshape(Ysol(:,iupmin),Nx,Ny); y = Y(:,end);
    plot(x,y,'b--','LineWidth', 1);
%     Y = reshape(Ysol(:,ilowmax),Nx,Ny); y = Y(:,1);
    Y = reshape(Ysol(:,iupmax),Nx,Ny); y = Y(:,end);
    plot(x,y,'b--','LineWidth', 1);  hold off
%     leg = [leg;strcat('t',num2str(k))];
% axis([2 6 min(y)-0.01 -0.5+0.01])
axis([0.5 5.5 0.5-0.01 max(y)+0.01])
% legend(leg)
% set(gca, 'fontsize', 15);
% plot(t(round([10:5:35]*Nt/50)),u0(round([10:5:35]*Nt/50)),'ok' ...
%    , 'MarkerFaceColor', 'k', 'MarkerSize', 10);
% text(t(round([10:5:35]*Nt/50)),u0(round([10:5:35]*Nt/50)), ...
%     ['(t_1)';'(t_2)';'(t_3)';'(t_4)';'(t_5)';'(t_6)'], 'HorizontalAlignment', 'right' ...
%     , 'VerticalAlignment', 'bottom', 'FontSize', 30);
% hold off; axis([5 37 min(u0)-0.1 max(u0)+0.1])
% % axis([10 50 min(u0)-0.01 max(u0)+0.01])
% set(gca, 'fontsize', 15);

h=gcf;
% % set(h,'renderer','opengl');
set(h, 'PaperPosition', [0 0 10 6])
% % [0 0 10 6],[0 0 30 24],[-1.5 0 23 18]
set(h,'PaperSize', [10 6]);
% % [10 6],[30 24],[19.5 18]

% print(h,'PulsNSymTwoElaAneurysmu0','-dpdf')
% print(h,'PulsNSymTwoElaAneurysmT6SL','-dpdf')
% print(h,'PulsNSymTwoElaAneurysmT1P','-dpdf')
print(h,'WallDefUP','-dpdf')
% print(h,'WallDefLOW','-dpdf')