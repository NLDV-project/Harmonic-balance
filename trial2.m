function trial2()
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
X=(1-wi.^2).^-1;
xw = fft(F);
x=fminsearch(logic(),xw);
global wi, f, n, T, i,wi, F,X,xw

function res = logic()
aw = wi.^2.*xw;
vw = wi.*xw;
at = ifft(aw);
vt = ifft(vw);
xt = at./wi.^2;
res = @(xw) sum(abs(at+xt+0*vt-F)); 
end

figure(1)
S = F./(1-w^2);
s1 = ifft(x)/(1-w^2);
plot(t,s1,t,S,'*');
end