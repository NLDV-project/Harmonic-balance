%Program for learning the Fast Fourier Transforms application on a simple harmonic balance problem x"+x = sin(wd*t), wd = Driving Frequency
function x = NLDV()
clc
clear all
close all
%Assumed frequency
f =3;

% initial guess of displacement for the equation
prompt = 'your initial or start position angle';
A= input(prompt);
prompt1 = 'Total period of single continuous cycle';
p=(input(prompt1))/4;
%Displacement
x = [sin(A+(p*0)) sin(A+(p*1)) sin(A+(p*2))]
%converting theta from degrees to radians
theta = [p*0 pi*((p*1)/180) pi*((p*2)/180)]
%Applying fast fourier transformation for finding closest Displacement value
y1 = fft(x)

% time
T = inv(f);
delta_t = T/4
t = delta_t*[0,1,2]
i = sqrt(-1);
w_i = (pi-(theta.*delta_t/(2*pi)))*i
% w_i = [0,3*i,6*i,-3*i];
v1 = w_i.*y1;
u1 = ifft(v1);

a = init();
for n=1:1:3
        x2(n) = 2*sin(2*pi*n*t(1,n))-x(n);
    end
    A = fft(x2);
    X(1) = 0;
    c = w_i.^2;
    for n=2:1:3
        X(n) = (-1/(w_i(n)^2))*A(n);
    end
x = fft(X);
d1 = x;
count = 1;
while d1-a~= 0
    a= d1;
    for n=1:1:3
        x2(n) = 2*sin(2*pi*n*t(1,n))-x(n);
    end
    A = fft(x2);
    X(1) = 0;
    c = w_i.^2;
    for n=2:1:3
        X(n) = (-1/(w_i(n)^2))*A(n);
    end
x = fft(X);
d1 = x;
count = count+1
end 
end 

function x = init()
f = 3;
% initial guess
x = [sin(0) sin(120) sin(240)]
theta = [0 120/180*pi 240/180*pi]
y1 = fft(x);

% time
T = inv(f);
delta_t = T/4;
t = delta_t*[0,1,2];
i = sqrt(-1);
w = 2*pi*f;
w_i = (pi-(theta.*delta_t/(2*pi)))*i
v1 = w_i.*y1;
u1 = ifft(v1);
for n=1:1:3
x2(n) = sin(2*pi*n*t(1,n))-x(n);
end
A = fft(x2);
X(1) = 0;
c = w_i.^2;
for n=2:1:3
    X(n) = (-1/(w_i(n)^2))*A(n);
end
x = fft(X);
end


% function x = ite()
% f=3;
% T = inv(f);
% delta_t = T/4;
% t = delta_t*[0,1,2,3];
% for n=1:1:4
% x2(n) = sin(2*pi*n*t(1,n))-x(n);
% end
% A = fft(x2);
% X(1) = 0;
% c = w_i.^2;
% for n=2:1:4
%     X(n) = (-1/(w_i(n)^2))*A(n);
% end
% x = fft(X);
% end

% function [d x] = ite2()
% d = d1(NLDV) - ite()
% end


