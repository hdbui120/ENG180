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

[a1,a2,a3,a4] = linsplin(del_1,k,n_1);
[b1,b2,b3,b4] = linsplin(del_2,k,n_2);
[c1,c2,c3,c4] = linsplin(del_3,k,n_3);

figure(1);
subplot(2,2,1);
hold on; grid on;
linsplinplot(a1,a2,a3,a4);
subplot(2,2,2);
hold on; grid on;
linsplinplot(b1,b2,b3,b4);
subplot(2,2,[3 4]);
hold on; grid on;
linsplinplot(c1,c2,c3,c4);

%plotting function
function linsplinplot(a,b,c,n)
    %plotting
    plot(b, c, 'bo',b, c, '--')
    for i = 1:n-1
        plot(a{i}(1,:),a{i}(2,:),'.k')
    end
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
        theta(i+1) = theta(i) + del*(k^i);
    end

    %generating range points (change orginal eqn here)
    f = sin(theta./4).^3;
    
    %setting up 10 points for each spline
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
