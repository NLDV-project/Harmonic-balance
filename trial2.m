function trial2()
clc
close all
clear all
%x''+x = sin(0.9t)
w=2;
f = w/(2*pi);
n=31;
T=inv(f);
t=linspace(0,T,n)
i = sqrt(-1);
wi= [0,(1:1:((n-1)/2)),-((n-1)/2:-1:1)]*((2*pi/T))*i;
F= cos(w.*t) % forcing function
f_n = 0.1*(cos(w*t)).^3+(0.2*w*sin(w.*t)) % non linear force
f = (F-f_n)/(1-w^2)+eps;
xw = fft(f);
x=fminsearch(logic(),xw);
global wi, f, n, T, i,wi, F,xw;

function res = logic()
aw = wi.^2.*xw;
vw = wi.*xw;
at = ifft(aw);
vt = ifft(vw);
xt = at./wi.^2;
res = @(xw)(sum(abs((at+xt-F+f_n)))); 
end

% analytical solution
%x¨(t) + ?x?(t) + ?x(t) + ?x(t)^3  = ? cos(?t)
function ydiff=duffing(t,y)
alpha=0.1;
beta=1;
gamma=1;
delta=0.2;
dydt=[y(2); -delta*y(2)-beta*y(1)-alpha*(y(1)^3)+gamma*cos(w*t)];
ydiff=dydt;
end

function ydiff=duff(t,y)
alpha=0.1;
beta=1;
gamma=1;
delta=0.2;
dydt=[y(2); -delta*y(2)-beta*y(1)-alpha*(y(1)^3)];
ydiff=dydt;
end

figure(1)
s1 = ifft(x)
plot(t,s1)

figure(2)
[t,y2]=ode45(@duffing,t,[0;0])
[t,y1]=ode45(@duff,t,[0;0])
y = y2-y1
plot(t,y,t,s1);

end