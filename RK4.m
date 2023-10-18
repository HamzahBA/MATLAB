function [t,y] = RK4(f,t,y0)

y = [y0'];
for i=1:length(t)-1
    h = t(i+1)-t(i);
    k1 = f(t(i),y0);
    k2 = f(t(i)+h/2,y0+h/2*k1);
    k3 = f(t(i)+h/2,y0+h/2*k2);
    k4 = f(t(i)+h,y0+h*k3);
    
    y0 = y0 + h/6*(k1+2*k2+2*k3+k4);
    y=[y;y0'];
end
end