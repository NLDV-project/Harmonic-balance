clc
close all
clear all
%x''+x = sin(0.9t)
w=0.9;
f = w/(2*pi);
n=31;
T=inv(f);
delta_t = T/(n-1);
t=linspace(0,T,n);
i = sqrt(-1);
wi= [0,(1:1:((n-1)/2)),-((n-1)/2:-1:1)]*((2*pi/T))*i;
F= sin(w.*t);
X=F./(1-w.^2);

x_w= fft(F);
vw = wi.*x_w; 
vt = ifft(vw);
aw = (F-x_w);
xw = aw./wi.^2;
xw(1) = 0;
xt = ifft(xw);
at = wi.^2.*xt;
vt = wi.*xt;
for n1 =1:1:(n-1)/2
    res(n1) =abs(at(n1)+F(n1)-xt(1)-((0.5*xt(n1+1))*exp(i*w*t(n1)))-((0.5*xt(end-n1))*exp(i*w*t(n1))));
end
e = sum(res);
x = xt./(1-w.^2);
figure(1)
plot(t,x,t,X,'*');



