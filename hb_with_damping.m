%equation x"+x'+x=sinwt
function hb_with_damping()
clc
clear all
close all
global N T w t F W X
N=11; %No. of points assumed
T=65; %Total time
w=2*pi/T;
t=linspace(0,T,N+1);%actually they will be N equally divided points in continuous system
t=t(1:end-1);%we required N-1 points in the current system
F=sin(w*t);%forcing term
iw=(0:ceil(N-1)/2)*1i*(w); miw=(-1i)*(floor(N/2):-1:1)*(w); %- & + values of w
W=[iw,miw];
X=fft(F);

x=fminsearch(error(), X);

function residue=error()
dotX=W.*X;
dotx=ifft(dotX);
ddotX=(W.^2).*X;
ddotx=ifft(ddotX);
x0=ddotx/W.^2; %displacement in time domain

residue=@(X) sum(abs(ddotx+dotx+x0-F).^2)

end
 figure(1)
 a=w*sqrt(1-(ifft(X)).^2)-(ifft(X)).*(1-w^2);
 b=w*sqrt(1-F.^2)-F.*(1-w^2);
 plot(t,a./(w^2-w^4-1),t,b./(w^2-w^4-1),'*-')
 legend('x-fft','x-analytical')
 x_0=-0.1;v0=0;
 fnc = @(t,x)[x(2);sin(w*t)-x(2)-x(1)]
 [tspan,x_ode] = ode45(fnc,t,[x_0 v0])
 figure(2)
 plot(t,x_ode(:,1),t,a./(w^2-w^4-1),'o-')
 legend('x-ode','x-fft')
end
