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

% Set up matrices for cramer's
Al = [(x(2)-x(1)),((x(2)-x(1))^2);(x(3)-x(1)),((x(3)-x(1))^2)];
bl = [y(2);y(3)];
Ar = [(x(n-1)-x(n)),((x(n-1)-x(n))^2);(x(n-2)-x(n)),((x(n-2)-x(n))^2)];
br = [(y(n-1)-y(n));(y(n-2)-y(n))];
k1 = cramer(Al,bl,1);
kn = cramer(Ar,br,1);

% making tridiagonal for Thomas3 and solve for k1 to kn
a = (1/6)*ones(1,n);
b = (4/6)*ones(1,n);
c = (1/6)*ones(1,n);
e = ones(1,n);
for i = 2:n-1
    e(i) = (y(i+1)-y(i-1))/(2*(x(i+1)-x(i))); 
end
e(1) = k1;
e(n) = kn;
k = THOMAS3(a,b,c,e,n);
k(1) = k1;
k(n) = kn;


for i = 1:n-1
    m(i) = (3*(y(i+1)-y(i))/(x(i+1)-x(i))^2)-((k(i+1)+2*k(i))/(x(i+1)-x(i)));
    d(i) = ((k(i+1)+k(i))/(x(i+1)-x(i))^2)-((2*(y(i+1)-y(i)))/((x(i+1)-x(i))^3));
end

splines = cell(n-1,1);
for i = 1:n-1
    xspline = linspace(x(i),x(i+1),10);
    yspline = y(i)+k(i)*(xspline-x(i))+m(i)*(xspline-x(i)).^2+d(i)*((xspline-x(i)).^3);
    splines{i} = [xspline;yspline];
end

figure(1);
plot(x, y, 'bo')
hold on; grid on;
for i = 1:n-1
    plot(splines{i}(1,:),splines{i}(2,:),'.k')
end
title('Left Quadratic Spline Interpolation');
xlabel('Theta');
ylabel('f','Rotation',0);

%take matrix a,b, and column values then use cramer's rule 
function cram = cramer(a,b,col)
    deterA = (a(1,1).*a(2,2))-(a(1,2).*a(2,1));
    ai = a;
    ai(:,col) = b;
    deterAI = (ai(1,1).*ai(2,2))-(ai(1,2).*ai(2,1));
    cram = deterAI/deterA;
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