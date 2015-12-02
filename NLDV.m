function x = NLDV()
clc
clear all
close all

f =3;

% initial guess
x = [sin(0) sin(120) sin(240)]
theta = [0 120/180*pi 240/180*pi]
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


