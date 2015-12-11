%Linear forcing differential equation x"+x=sin(wt)
function harmonic_balance()
clc
clear all
close all
global N T w t F W X
N=5;
T=2;
w=2*pi/T;
t=linspace(0,T,N+1);%actually they will be N equally divided points in continuous system
t=t(1:end-1);%we required N-1 points in the current system
F=sin(w*t);%forcing term
iw=(0:ceil(N-1)/2)*1i*(w);
miw=(-1i)*(floor(N/2):-1:1)*(w);
W=[iw,miw];
X=fft(F);
x=fminsearch(error(), X);

function residue=error()
dotX=W.*X;
dotx=ifft(dotX);
ddotX=(W.^2).*X;
ddotx=ifft(ddotX);
x0=ddotx/W.^2;

residue=@(X) sum(abs(ddotx+x0-F))

end
 figure(1)
 plot(t,ifft(x)/(1-w^2),t,F./(1-w^2),'*-')
 legend('x-fft','x-analytical')
 x_0=0;v0=0;
 fnc = @(t,x)[x(2);sin(w*t)-x(1)]
 [tspan,xval] = ode45(fnc,t,[x_0 v0])
 figure(2)
 plot(t,xval(:,1),t,ifft(x)/(1-w^2),'o-')
 legend('x-val','x-fft')
end
