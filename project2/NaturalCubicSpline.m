clear, clc; 
close all;

clear, clc; 
close all;

    %---scope---%
% 1. solve ditriagonal for k
% 2. solve m and d
% 3. make spline

% make data points
n = 10;
x = linspace(0,2*pi,n);
y = sin(x./4).^3;

% exact case
xexact = linspace(0,2*pi(),100);
yexact = sin(xexact./4).^3;

s1 = naturalCubic(x,y,n);
naturalCubicPlot(x,y,xexact,yexact,s1,n)

function splines = naturalCubic(x,y,n)
    % making tridiagonal for Thomas3 and solve for k1 to kn
    a = (1/6)*ones(1,n);
    b = (4/6)*ones(1,n);
    c = (1/6)*ones(1,n);
    e = ones(1,n);
    for i = 2:n-1
        e(i) = (y(i+1)-2*y(i)+y(i-1))/(x(i+1)-x(i));
    end
    e(1) = 0;
    e(n) = 0;
    m = THOMAS3(a,b,c,e,n);
    m(1) = 0;
    m(n) = 0;
    
    splines = cell(n-1,1);
    for i = 1:n-1
        xspline = linspace(x(i),x(i+1),10);
        deltax = x(i+1)-x(i);
        a0 = m(i)/(6*deltax);
        a1 = m(i+1)/(6*deltax);
        a2 = y(i)/deltax-(m(i)*deltax)/6;
        a3 = y(i+1)/deltax-(m(i+1)*deltax)/6;
        yspline = a0.*(x(i+1)-xspline).^3 + a1.*(xspline-x(i)).^3 + a2.*(x(i+1)-xspline) + a3.*(xspline-x(i));
        splines{i} = [xspline;yspline];
    end
end

function naturalCubicPlot(x,y,xexact,yexact,splines,n)
    figure(1);
    plot(x, y, 'bo', xexact, yexact, '--')
    hold on; grid on;
    for i = 1:n-1
        plot(splines{i}(1,:),splines{i}(2,:),'.k')
    end
    title('Natural Cubic Spline Interpolation');
    xlabel('Theta');
    ylabel('f','Rotation',0);
end

function x = THOMAS3(a,b,c,d,n)

    %initial condition
    bbar(1) = b(1);
    cbar(1) = c(1);
    dbar(1) = d(1);

    %making upper triangle
    for i = 2:n
        multiplier = a(i)./bbar(i-1);
        abar(i) = a(i) - bbar(i-1).*multiplier;
        bbar(i) = b(i) - cbar(i-1).*multiplier;
        cbar(i) = c(i);
        dbar(i) = d(i) - dbar(i-1).*multiplier;
    end
    
    %initialize x of size n
    x = ones(1,n);

    %initialize end condition
    x(n) = dbar(n)/bbar(n);

    % Upward substitution AKA zip it up
    for i = n-1:-1:1
        x(i) = (dbar(i)-(cbar(i)*x(i+1)))/bbar(i);
    end
end