clc, clear; 
close all;

% Lord fogive me for what im about to do with this code
% make da grid
delta = .01; 
t = 0:delta:10;


[xref,vref,x,v] = ivp(1,0,delta,'eulerF',t);
[xref1,vref2,xbe,vbe] = ivp(1,0,delta,'eulerB',t);

% Plot area 
figure(1)
plot(t,x,t,xref)
hold on;
plot(t,xbe)
hold on;
xlabel('t')
ylabel('x')

figure(2)
plot(t,v,t,vref)
hold on;
plot(t,vbe)
hold on;
xlabel('t')
ylabel('v') 

function [xexact, vexact, x,v] = ivp(v0,x0,delta,method,t)
m = 3; g = 2;
c1 = m*g; c2 = m;
ydot = @(v) [c1/m-c2/m*v;v];
xexact = m/c2*(c1/c2-1)*exp(-c2/m*t)-m/c2*(c1/c2-1)+c1/c2*t;
vexact = -(c1/c2-1)*exp(-c2/m*t)+(c1/c2);

%Function handler for different methods
fEuler = @(dy,y,h) dy*h+y;
bEulerV = @(y,c1,c2,m,h) (y*m+c1*h)/(m+c2*h);  
bEulerX = @(y,v,h) y+v*h;

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
        gv = v0;
        gx = x0;
        v(1) = v0;
        x(1) = x0;
        % Fixed point iterations
        for i = 1:length(t)-1
            v(i+1) = bEulerV(gv,c1,c2,m,delta);
            x(i+1) = bEulerX(gx,v(i),delta);
            gv = v(i+1);
            gx = x(i+1);
        end

end

end
