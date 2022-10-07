clear, clc;
close all;

% 3 cases data
n1 = 10;
n2 = 20;
n3 = 40;
del1 = .4186;
del2 = .1106;
del3 = .0146;
k = 1.1;

%exact original fncs
xexact = linspace(0,2*pi(),100);
yexact = sin(xexact./4).^3;

%running methods with inputs from 3 cases
[s1,x1,y1] = LQuadSpline(del1,k,n1);
[s2,x2,y2] = LQuadSpline(del2,k,n2);
[s3,x3,y3] = LQuadSpline(del3,k,n3);

% Plotting
figure;
subplot(2,2,1);
LQuadPlot(x1,y1,xexact,yexact,s1,n1,'Case 1'); 
subplot(2,2,2);
LQuadPlot(x2,y2,xexact,yexact,s2,n2,'Case 2'); 
subplot(2,2,[3 4]);
LQuadPlot(x3,y3,xexact,yexact,s3,n3,'Case 3'); hold on;

%setting up 10 points for each spline
function [splines,x,y] = LQuadSpline(del,kstr,n)
   
    %generating grid
    x = ones(1,n);
    x(1) = 0;
    x(n) = 2*pi();
    for i = 1:n-2
        x(i+1) = x(i) + del*(kstr^(i-1));
    end
    
    %generating grid range points (change orginal eqn here)
    y = sin(x./4).^3;
    
    A = [(x(2)-x(1)),((x(2)-x(1))^2);(x(3)-x(1)),((x(3)-x(1))^2)];
    b = [y(2);y(3)];
    k = ones(n-1,1);
    m = ones(n-1,1);
    k(1) = cramer(A,b,1);
    
    for i = 1:n-1
        m(i) = ((y(i+1)-y(i))-k(i)*(x(i+1)-x(i)))/((x(i+1)-x(i))^2);
        k(i+1) = k(i) + 2*m(i)*(x(i+1)-x(i));
    end

    splines = cell(n-1,1);
    for i = 1:n-1
        xspline = linspace(x(i),x(i+1),10);
        yspline = y(i)+k(i)*(xspline-x(i))+m(i)*(xspline-x(i)).^2;
        splines{i} = [xspline;yspline];
    end
end

% plot x-y data points
% plot exact function
% plot splines
function LQuadPlot(x,y,xexact,yexact,s,n,str)
    plot(x, y, 'bo',xexact,yexact, '--')
    hold on; grid on;
    for i = 1:n-1
        plot(s{i}(1,:),s{i}(2,:),'.k')
    end
    title(str);
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

