function trial2()
clc
close all
clear all
%x''+x = sin(0.9t)
w=0.9;
f = w/(2*pi);
n=3;
T=inv(f);
t=linspace(0,T,n);
i = sqrt(-1);
wi= [0,(1:1:((n-1)/2)),-((n-1)/2:-1:1)]*((2*pi/T))*i;
F= sin(w.*t); % forcing function
f_n = 0;
f = (f_n+F)./(1-w.^2);
xw = fft(f);
x=fminsearch(logic(),xw);
global wi, f, n, T, i,wi, F,xw, f_n

function res = logic()
aw = wi.^2.*xw;
vw = wi.*xw;
at = ifft(aw);
vt = ifft(vw);
xt = at./wi.^2;
res = @(xw)(sum(abs((at+f+0*vt-F).^2))); 
end

figure(1)
S = F./(1-w.^2);
s1 = ifft(x);
plot(t,s1,t,S,'*');
end