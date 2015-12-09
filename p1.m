clc
clear all
close all

Fs=1;    % Sample Frequency
N=99;     % Sample Points
Ts=1/Fs; % Sample Time
dt=Ts/N   % Time sampling interval
t=(0:(N-1))*dt;
x=sin(2*pi*Fs*t);
X=(fft(x));

w_i=[0:(N-1)/2,((N-1)/2):-1:1]*(2*pi*sqrt(-1)*Fs);
X_1=X%./(1-(w_i.^2));
x_1=ifft(X_1);
plot(t,(x_1),t,x,'*')
xlabel('time')
ylabel('Amplitude')


