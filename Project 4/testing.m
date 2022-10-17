clc, clear; 
close all;

k = 8;

c1 = 1.6;
c2 = 1.1;
%     a = -.13;
%     b = .44;
%     c = .25;
%     alpha = 2;
x = ones(k,1);
r = ones(k,1);
x(1) = 1.6;

for i=1:k-1
    r(i) = (c2/c1*x(i))^(1/2.1);
    rdot = (c2/c1/2.1)*(c2/c1*x(i))^((1-2.1)/2.1);
    x(i+1) = x(i) - r(i)/rdot;
end
r(k) = (c2/c1*x(k))^(1/2.1);

y = linspace(-2,2,100);
rexact = (1.6/1.1).*y.^2;
figure(1);
plot(y,y,y,rexact)
xline(x(1))
hold on; grid on;
plot(x(k),r(k-1),'bo')
hold on;