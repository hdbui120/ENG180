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
figure('Name','Fixed Point Iterations','NumberTitle','off')
subplot(3,1,1)
[r11,r21]=p2f('a',k,.5,3,.7);
title(sprintf('Propellant Burnt when n = .5, r1 = %.2f, r2 = %.2f', r11, r21))
xlabel('x')
ylabel('y')
hold on;
subplot(3,1,2)
[r12,r22]=p2f('a',k,2.1,1,.6);
title(sprintf('Propellant Burnt when n = 2.1, r1 = %.2f, r2 = %.2f', r12,r22))
xlabel('x')
ylabel('y')
hold on;
subplot(3,1,3)
[r13,r23]=p2f('b',k,1,2,.8);
title(sprintf('Semenov Equation, r1 = %.2f, r2 = %.2f', r13, r23))
xlabel('x')
ylabel('y')
hold on;

% figure('Name','Secant Method',NumberTitle='off')
% [root1,root2] = p2secant('a',1,1.3,.2,.2,30,.5);
% title(sprintf('Propellant Burnt when n = .5, r1 = %.2f, r2 = %.2f', root1, root2))
% xlabel('x')
% ylabel('y')

% subplot(4,1,4)
% [newtonR11, newtonR22]=p2n('a',k,.5,3,2);
% title('Propellant Burnt when n=.5 using Newton')
% xlabel('x')
% ylabel('y')
% hold on;


% Fixed point iterations problem 2A
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
    x(1) = g;
    if strcmp(e,'a') 
        %===========root 1===========%
        fun1 = @(x)(c1/c2).*x.^n; 
        fun2 = @(x)(c2/c1).*x.^(1/n);
        for i = 1:k-1
            r(i) =  fun1(x(i));%propellant burning root equation 1
            x(i+1) = r(i);
        end
        r(k) = fun1(x(k));
        y = linspace(0,s,100);
        rexact = fun1(y);
    
        plot(y,y,y,rexact,'--')
        hold on;
        xline(g);
        plot(x(k),r(k),'^','MarkerFaceColor','red')
        hold on;
        r1 = r(k);
        %===========root 2===========%
        for i = 1:k-1
            r(i) = fun2(x(i)); %propellant burning root equation 2
            x(i+1) = r(i);
        end
        r(k) = fun2(x(k));
        y2 = linspace(0,s,100);
        rexact2 = fun2(y);
    
        plot(y2,rexact2,'--')
        hold on;
        plot(x(k),r(k),'^','MarkerFaceColor','red')
        hold on;
        r2 = r(k);        
    elseif strcmp(e,'b')
        %semenov root equation 1
        semenov1 = @(x) (c.*exp(1).^(-alpha.*x)-a)./b;
        for i = 1:k
            r(i) = semenov1(x(i));
            x(i+1) = r(i);
        end
        r(k) = semenov1(x(k));
        y = linspace(0,s,20);
        rexact = semenov1(y);
    
        plot(y,y,y,rexact,'--')
        xline(g);
        hold on;
        plot(x(k),r(k),'^','MarkerFaceColor','red')
        hold on;
        r1 = r(k);
        
        %semenov root equation 2
        semenov2 = @(x) alpha./(log((alpha+b.*x)./c));
        for i = 1:k
            r(i) =  semenov2(x(i));
            x(i+1) = r(i);
            if semenov2(x(i+1))-r(i)>1
                disp("pass tolerance")
                break
            end
        end
        r(k) = semenov2(x(k));
        y = linspace(0,s,20);
        rexact2 = semenov2(y);

        plot(y,y,y,rexact2,'--')
        xline(g);
        hold on;
        plot(x(k),r(k),'^','MarkerFaceColor','red')
        hold on;
        r2 = r(k);
    end
end

%Approximate with newton method
function [r1,r2] = p2n(e,k,n,s,g)

c1 = 1.6;
c2 = 1.1;
a = -.13;
b = .44;
c = .25;
alpha = 2;

%root function 1 for propellant burnt
if strcmp(e,'a')
    x = ones(k,1);
    r = ones(k,1);
    x(1) = g;
    for i=1:k-1
        r(i) = (c1/c2)*x(i)^n; 
        rdot = (n*c1/c2)*x(i)^(n-1);
        x(i+1) = x(i) - r(i)/rdot;
    end
    r(k) = (c1/c2)*x(k)^n;
    
    y = linspace(-2,s,100);
    rexact = (c1/c2).*y.^n;
    
    plot(y,rexact,'--')
    xline(x(1))
    hold on; grid on;
    plot(x(k),r(k),'bo')
    hold on;
    r1 = r(k);

    %root function 2 for propellant burnt
    for i=1:k-1
        r(i) = (c2/c1*x(i))^(1/n);
        rdot = (c2/c1/n)*(c2/c1*x(i))^((1-n)/n);
        x(i+1) = x(i) - r(i)/rdot;
    end
    r(k) = (c2/c1*x(k))^(1/n);
    
    y = linspace(-2,s,100);
    rexact = (c2/c1).*y.^(1/n);
    
    plot(y,rexact,'--')
    xline(x(1))
    hold on; grid on;
    plot(x(k),r(k),'bo')
    hold on;
    r2 = r(k);
elseif strcmp(e,'b')
    x = ones(k,1);
    r = ones(k,1);
    x(1) = g;
    %====root equation 1 for semenov====%
    for i = 1:k-1
        r(i) = (c*exp(1)^(-alpha/x(i))-a)/b;
        rdot = alpha*c/b*x(i)^2*e^(-alpha/x(i));
        x(i+1) = x(i) - r(i)/rdot;
    end
    r(k) = (c*exp(1)^(-alpha/x(k))-a)/b;
    
    y = linspace(-2,s,100);
    rexact = (c.*exp(1).^(-alpha./y)-a)./b;
    
    plot(y,rexact,'--')
    xline(g)
    hold on; grid on;
    plot(x(k),r(k),'bo')
    hold on;
    r1 = r(k);
    
    %====root 2 for semenov====%
    for i = 1:k-1
        r(i) = alpha*(log((alpha+b*x(i))/c))^-1;
        rdot = alpha*b/c*(log((alpha+b*x(i))/c))^-2;
        x(i+1) = x(i) - r(i)/rdot;
    end
    r(k) = alpha*log((alpha+b*x(k))/c)^-1;
    y = linspace(-2,s,100);
    rexact2 = -alpha.*log((alpha+b.*y)./c).^-1;
    
    plot(y,rexact2,'--')
    xline(g)
    hold on; grid on;
    plot(x(k),r(k),'bo')
    hold on;
    r2 = r(k);
end
end

function [root1,root2] = p2secant(eqn,x1,x2,x11,x22,k,n)
c1 = 1.6;
c2 = 1.1;
a = -.13;
b = .44;
c = .25;
alpha = 2;

if strcmp(eqn, 'a')
    rootFunc1 = @(x) (c1/c2)*x.^n;
    rootFunc2 = @(x) ((c2/c1)*x).^(1/n);
    x(1) = x1;
    x(2) = x2;
    
    root = ones(k,1);
    root(1) = rootFunc1(x(1));
    
    for i = 2:k-1
        root(i) = rootFunc1(x(i));
        rdot = root(i)-root(i-1)/(x(i)-x(i-1));
        x(i+1) = x(i) -root(i)/rdot;
    end
    root(k)=rootFunc1(x(k));
    root1 = root(k);
    y = linspace(0,5,20);
    xline(x(1))
    hold on; grid on;
    plot(y,rootFunc1(y),x(k),root(k),'bo')
    
    x(1) = x11;
    x(2) = x22;
    root(1) = rootFunc2(x(1));
    for i = 2:k-1
        root(i) = rootFunc2(x(i));
        rdot = root(i)-root(i-1)/(x(i)-x(i-1));
        x(i+1) = x(i) -root(i)/rdot;
    end
    root(k) = rootFunc2(x(k));
    root2 = root(k);
    y = linspace(0,5,20);
    xline(x(1))
    hold on; grid on;
    plot(y,rootFunc2(y),x(k),root(k),'bo')

elseif strcmp(eqn,'b')
end

end