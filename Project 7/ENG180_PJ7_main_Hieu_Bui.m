clc, clear; 
close all;

% Lord fogive me for what im about to do with this code
% make da grid
delta = .01; 
t = 0:delta:10;


[xref,vref,x,v] = ivp(1,0,delta,'eulerF',t);

% Plot area 
figure(1)
plot(t,x,t(1:30:length(t)),xref(1:30:length(t)),'ko')
xlabel('t')
ylabel('x')
hold on; grid on;
figure(2)
plot(t,v,t(1:30:length(t)),vref(1:30:length(t)),'ko')
xlabel('t')
ylabel('v')
hold on; grid on;

function [xexact, vexact, x,v] = ivp(v0,x0,delta,method,t)
m = 3; g = 2;
c1 = m*g; c2 = m;
ydot = @(v) [c1/m-c2/m*v;v];
fEuler = @(dy,y,h) dy*h+y;
xexact = m/c2*(c1/c2-1)*exp(-c2/m*t)-m/c2*(c1/c2-1)+c1/c2*t;
vexact = -(c1/c2-1)*exp(-c2/m*t)+(c1/c2);

switch method
    case 'eulerF'
        yold = [v0;x0];
        for i = 1:length(t)-1
            y = yold(:,i);
            dy = ydot(y(1,1));
            ynew = fEuler(dy,y,delta);
            yold(:,i+1) = ynew;
        end
        y = yold;
        v = y(1,:);
        x = y(2,:);
    case 'eulerB'
end

end
