%dotdot x+x=cos(2*w*t)
clc
clear all
close all
f=3;
w=2*pi*f;
T=1/f;
N=100;
t=(0:1:N-1)*T/N;
if (mod(N,2))==0
theta= ((2*pi*2/((N)))*linspace(ceil(-N/2),ceil(N/2),1));
else
    theta= ((2*pi*2/((N+1)))*linspace(ceil(-N/2),ceil(N/2),1));
end
x = cos(2*w*t); %initial guess (time domain)
X=(fft(x)) %frequency domain
w1 = theta/T;
X2=X./(1-(4*w1.^2));
X1=(ifft(X2))
x_theo= (1/(1-4*w^2))*cos(2*w*t)
figure (1)
plot(t,X1,'r')
hold on
plot(t,x_theo,'--')
