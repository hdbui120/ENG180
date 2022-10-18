clc, clear; 
close all;

k = 8;

c1 = 1.6;
c2 = 1.1;
a = -.13;
b = .44;
c = .25;
alpha = 2;
x = ones(k,1);
r = ones(k,1);
x(1) = .2;
n = 2;

for i = 1:k-1
    r(i) = (c*exp(1)^(-alpha/x(i))-a)/b;
    rdot = (alpha*c/(b*x(i)^2))*exp(1)^(-alpha/x(i));
    x(i+1) = x(i) - r(i)/rdot;
end
r(k) = (c*exp(1)^(-alpha/x(k))-a)/b;

y = linspace(-2,2,100);
rexact = (c.*exp(1).^(-alpha./y)-a)./b;

figure(1);
plot(y,y,y,rexact)
xline(x(1))
hold on; grid on;
% plot(x(k),r(k-1),'bo')
% hold on;