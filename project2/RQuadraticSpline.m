clear, clc;
close all;

%case 1
n = 10;
del = .4186;
k = 1.1;

%generating grid
x = ones(1,n);
x(1) = 0;
x(n) = 2*pi();
for i = 1:n-2
    x(i+1) = x(i) + del*(k^(i-1));
end

%generating grid range points (change orginal eqn here)
y = sin(x./4).^3;

%exact original fncs
xexact = linspace(0,2*pi(),100);
yexact = sin(xexact./4).^3;

%running methods with inputs from case 1
case1 = RQuadSpline(x,y,n);
RQuadPlot(x,y,xexact,yexact,case1,n);



%setting up 10 points for each spline
function splines = RQuadSpline(x,y,n)   
    A = [(x(n-1)-x(n)),((x(n-1)-x(n))^2);(x(n-2)-x(n)),((x(n-2)-x(n))^2)];
    b = [(y(n-1)-y(n));(y(n-2)-y(n))];
    k = ones(n-1,1);
    m = ones(n-1,1);
    k(n) = cramer(A,b,1);
    
    for i = n:-1:2
        m(i) = ((y(i-1)-y(i))-k(i)*(x(i-1)-x(i)))/((x(i-1)-x(i))^2);
        k(i-1) = k(i) + 2*m(i)*(x(i-1)-x(i));
    end

    splines = cell(n-1,1);
    for i = 2:n
        xspline = linspace(x(i-1),x(i),10);
        yspline = y(i)+k(i)*(xspline-x(i))+m(i)*(xspline-x(i)).^2;
        splines{i} = [xspline;yspline];
    end
end

% plot x-y data points
% plot exact function
% plot splines
function RQuadPlot(x,y,xexact,yexact,s,n)
    figure(1);
    plot(x, y, 'bo',xexact,yexact, '--')
    hold on; grid on;
    for i = 2:n
        plot(s{i}(1,:),s{i}(2,:),'.k')
    end
    title('Quadratic Spline Interpolation');
    xlabel('Theta');
    ylabel('f','Rotation',0);
end

%take matrix a,b, and column values then use cramer's rule 
function cram = cramer(a,b,col)
    deterA = (a(1,1).*a(2,2))-(a(1,2).*a(2,1));
    ai = a;
    ai(:,col) = b;
    deterAI = (ai(1,1).*ai(2,2))-(ai(1,2).*ai(2,1));
    cram = deterAI/deterA;
end

