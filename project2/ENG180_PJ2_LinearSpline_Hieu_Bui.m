clear, clc; 
close all; 

%constant variables
n_1 = 10;
del_1 = .4186;
n_2 = 20;
del_2 = .1106;
n_3 = 40;
del_3 = .0146;
k = 1.1;

%exact original fncs
xexact = linspace(0,2*pi(),100);
yexact = sin(xexact./4).^3;

[a1,a2,a3,a4] = linsplin(del_1,k,n_1); %case 1
[b1,b2,b3,b4] = linsplin(del_2,k,n_2); %case 2
[c1,c2,c3,c4] = linsplin(del_3,k,n_3); %case 3

%integrate
area1 = splineArea(a1,n_1); 
area2 = splineArea(b1,n_2);
area3 = splineArea(c1,n_3);

%differentiate
diff1 = deriv(a1);
diff2 = deriv(b1);
diff3 = deriv(c1);

%plotting all 3 cases
figure(1);
subplot(2,2,1);
hold on; grid on;
linsplinplot(xexact,yexact,a1,a2,a3,a4,'Case 1');
subplot(2,2,2);
hold on; grid on;
linsplinplot(xexact,yexact,b1,b2,b3,b4,'Case 2');
subplot(2,2,[3 4]);
hold on; grid on;
linsplinplot(xexact,yexact,c1,c2,c3,c4,'Case 3');

%plotting function
function linsplinplot(xexact,yexact,a,b,c,n,str)
    %plotting
    plot(b, c, 'bo', xexact, yexact, '--')
    for i = 1:n-1
        plot(a{i}(1,:),a{i}(2,:),'.k')
    end
    title(str);
    xlabel('Theta');
    ylabel('f','Rotation',0);
end


%linear spline function. 
%The function returns splines, each with 10 discrete points per spline.
%Return data points and number of points
function [a,b,c,d] = linsplin(del, k, n)
    %generating domain points
    theta = ones(1,n);
    theta(1) = 0;
    theta(n) = 2*pi();
    
    for i = 1:n-2
        theta(i+1) = theta(i) + del*(k^(i-1));
    end

    %generating range points (change orginal eqn here)
    f = sin(theta./4).^3;
    
    %setting up 10 points for each spline
    %splines are stored as cell array, which contains 10 points between
    %data points
    s = cell(n-1,1);
    for j = 2:n
        m = (f(j)-f(j-1))/(theta(j)-theta(j-1));
        x = linspace(theta(j-1),theta(j),10);
        b = f(j)-m*theta(j);
        y = m.*x + b;
        s{j-1} = [x;y];
    end
    a = s;
    b = theta;
    c = f;
    d = n;
end

%integrating function
function area = splineArea(splines,n)
    area = 0;
    for j=1:n-1
        for i = 1:9
            littleArea = splines{j}(2,i)*(splines{j}(1,i+1)-splines{j}(1,i));
        end
        area = area + littleArea;
    end
end

%differentiating function
function diff = deriv(splines)
    diff = ones(3,1);
    for i=1:3
        diff(i) = (splines{i}(2,i+1)-splines{i}(2,i))/(splines{i}(1,i+1)-splines{i}(1,i));
    end
end