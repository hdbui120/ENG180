clc, clear; 
close all;

%general strategy for secant method
% r = root equation;
% x(1) = initial guess;
% rdot = [r(x(i))-r(x(i-1))]/(x(i)-x(i-1))
% 
% x(i+1) = x(i) - r(i)/rdot for 2<i<k
% r(k) = (r(x(k)))

k = 30;
c1 = 1.6; c2 = 1.1; n = .5;

[r1,r2] = p2secant(1,1.1,.3,.3,k,n);

function [root1,root2] = p2secant(x11,x21,x12,x22,k,n)
c1 = 1.6;
c2 = 1.1;
a = -.13;
b = .44;
c = .25;
alpha = 2;

rootFunc1 = @(x) (c1/c2)*x.^n;
x(1) = x11;
x(2) = x21;

root = ones(k,1);
root(1) = rootFunc1(x(1));

for i = 2:k-1
    root(i) = rootFunc1(x(i));
    rdot = root(i)-root(i-1)/(x(i)-x(i-1));
    x(i+1) = x(i) -root(i)/rdot;
end
root(k)=rootFunc1(x(k));
root = root(k);

y = linspace(0,5,20);
figure(1)
xline(x(1))
hold on;
plot(y,rootFunc1(y),x(k),root(k),'bo')

end