%{
    ENG 180 Project 4
    File:           ENG180_PJ3_main_Hieu_Bui.m
    Author:         Hieu Bui
    Date:           10/19/2022
    Description:    Finding answers through iterative approximation

%}

clc, close all;
clear;

k = 8;
figure(2)
subplot(3,1,1)
r2a1=p2f('a',k,.5,3,.7);
title('Propellant Burnt when n = .5')
xlabel('x')
ylabel('y')
hold on;
subplot(3,1,2)
r2a2=p2f('a',k,2.1,1,.6);
title('Propellant Burnt when n = 2.1')
xlabel('x')
ylabel('y')
hold on;
subplot(3,1,3)
r2b=p2f('b',k,1,.7,.6);
title('Semenov Equation')
xlabel('x')
ylabel('y')
hold on;

% Main function for problem 2A
% e = which equation ['a' = 2a; 'b' = 2b]
% k is number of iterations to do
% n is the power for problem 2a
% s is number of points for exact graph
% g is initial guess
% r is the approximated root
function [r1,r2] = p2f(e,k,n,s,g)
    c1 = 1.6;
    c2 = 1.1;
    a = -.13;
    b = .44;
    c = .25;
    alpha = 2;
    r = zeros(k,1);
    x = ones(k,1);
    r(1) = g;
    x(1) = r(1);
    if strcmp(e,'a') 
        %===========root 1===========%
        for i = 1:k-1
            r(i) = (c1/c2)*x(i)^n; %propellant burning root equation 1
            x(i+1) = r(i);
        end
        r(k) = (c1/c2)*x(k)^n;
        y = linspace(0,s,100);
        rexact = (c1/c2).*y.^n;
    
        plot(y,y,y,rexact,'--')
        hold on;
        xline(g);
        plot(x(k),r(k),'^','MarkerFaceColor','red')
        hold on;
        for j = 2:length(r)
            plot([x(j-1) x(j-1)],[x(j-1) r(j-1)],'b--');
            pause(.1);
            plot([x(j-1) r(j-1)],[x(j) x(j)],'b--')
            pause(.1);
            hold on;
        end
        r1 = r(k);
        %===========root 2===========%
        for i = 1:k-1
            r(i) = (c2/c1)*x(i)^(1/n); %propellant burning root equation 2
            x(i+1) = r(i);
        end
        r(k) = (c2/c1)*x(k)^(1/n);
        y2 = linspace(0,s,100);
        rexact2 = (c2/c1).*y.^(1/n);
    
        plot(y2,rexact2,'--')
        hold on;
        plot(x(k),r(k),'^','MarkerFaceColor','red')
        hold on;
        for j = 2:length(r)
            plot([x(j-1) x(j-1)],[x(j-1) r(j-1)],'b--');
            pause(.1);
            plot([x(j-1) r(j-1)],[x(j) x(j)],'b--')
            pause(.1);
            hold on;
        end
        r2 = r(k);        
    elseif strcmp(e,'b')
        for i = 1:k
            r(i) = (c*exp(1)^(-alpha*x(i))-a)/b; %semenov root equation
            x(i+1) = r(i);
        end
        r(k) = (c*exp(1)^(-alpha*x(k))-a)/b;
        y = linspace(.4,s,100);
        rexact = (c.*exp(1).^(-alpha.*y)-a)./b;
    
        plot(y,y,y,rexact)
        xline(g);
        hold on;
        plot(x(k),r(k),'^','MarkerFaceColor','red')
        hold on;
        for j = 2:length(r)
            plot([x(j-1) x(j-1)],[x(j-1) r(j-1)],'b--');
            pause(.1);
            plot([x(j-1) r(j-1)],[x(j) x(j)],'b--')
            pause(.1);
            hold on;
        end
        r1 = r(k);
    end
end

function [r1,r2] = p2n(e,k,n,s,g)
    c1 = 1.6;
    c2 = 1.1;
    %     a = -.13;
    %     b = .44;
    %     c = .25;
    %     alpha = 2;

    x = ones(k,1);
    r = ones(k,1);

    for i=1:k-1
        r(i) = (c1/c2)*x(i)^n; %root function 1
        rdot = (n*c1/c2)*x(i)^(n-1);
        x(i+1) = x(i) - r(i)/rdot;
    end
    r(k) = (c1/c2)*x(k)^n;
    
    y = linspace(-2,2,100);
    rexact = (1.6/1.1).*y.^n;
    
    plot(y,y,y,rexact)
    xline(x(1))
    hold on; grid on;
    plot(x(k),r(k),'bo')
    hold on;
    r1 = r(k);

    for i=1:k-1
        r(i) = (c2/c1*x(i))^(1/n);
        rdot = (c2/c1/n)*(c2/c1*x(i))^((1-n)/n);
        x(i+1) = x(i) - r(i)/rdot;
    end
    r(k) = (c2/c1*x(k))^(1/n);
    
    y = linspace(-2,2,100);
    rexact = (c2/c1).*y.^(1/n);
    
    plot(y,y,y,rexact)
    xline(x(1))
    hold on; grid on;
    plot(x(k),r(k),'bo')
    hold on;
    r2 = r(k);
end