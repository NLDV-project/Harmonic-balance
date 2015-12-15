function nonlinear_damp()
clc
close all
clear all
w=0.08;
f = w/(2*pi);
n=15;
T=inv(f);
t=linspace(0,T,n)
i = sqrt(-1);
wi= [0,(1:1:((n-1)/2)),-((n-1)/2:-1:1)]*((2*pi/T))*i;
f = (cos(w.*t)+sin(w.*t))-(-w.*sin(w*t))/(1-w^2)+eps;
xw = fft(f);
x=fminsearch(logic(),xw);
global wi, f, n, T, i,xw;

function res = logic()
aw = wi.^2.*xw;
vw = wi.*xw;
at = ifft(aw);
vt = ifft(vw);
xt = at./wi.^2;
res = @(xw)(sum(abs((at+xt-(f*(1-w.^2)))))); 
end
% analytical solution
function ydiff=damping(t,y)
dydt=[y(2);-y(1)-w*sin(w*t)+(cos(w*t)+sin(w*t))];
ydiff=dydt;
end
function ydiff=damp(t,y)
dydt=[y(2); -y(1)-w*sin(w*t)];
ydiff=dydt;
end
s1 = ifft(x)
figure(1)
[t,y2]=ode45(@damping,t,[1;0])
[t,y1]=ode45(@damp,t,[0;0])
 y = y2-y1
plot(t,y(:,1),'*-',t,s1);
legend('x-ode','x-fft')
end