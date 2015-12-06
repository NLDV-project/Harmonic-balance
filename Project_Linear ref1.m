%dotdot x+x=sin(2*w*t)
clc
clear all
close all
f=3; w=2*pi*f; c=3; T=1/f; N=19;
linewidth=5
t=(0:1:N-1)*T/N;
theta= ((2*pi*2/((ceil(N))))*linspace(ceil(-N/2),ceil(N/2),1));
x = sin(c*w*t); %initial guess (time domain)
X=(fft(x)) %frequency domain
w1 = theta/T;
X2=X./(1-(c^2*w1.^2));
X1=(ifft(X2))
x_theo= (1/(1-c^2*w^2))*sin(c*w*t)
plot(t,X1,'r',t,x_theo,'g*--')