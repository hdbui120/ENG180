clear, clc;
close all;

%declaring and initializing variables
n = 4;
a = 2.*ones(1,n); % a(1) = 1
b = -5.*ones(1,n);
c = 2.*ones(1,n);
d = [-3,-1,-1,-3];

z = decomp(a,b,c,d,n)

function x = decomp(a,b,c,d,n)
%     initial condition
    beta(1) = b(1);
    gamma(1) = c(1);

    for i = 2:n
        alpha(i) = a(i)/beta(i-1);
        beta(i) = b(i) - alpha(i)*gamma(i-1);
        gamma(i) = c(i);
    end

%     initial condition for zip down
    y(1) = d(1);
%     zip down
    for i = 2:n
        y(i) = d(i) - alpha(i)*y(i-1);
    end
%     preallocate x
    x = ones(1,n);
%     initial condition for zipping up
    x(n) = y(n)/beta(n);
%     zip up
    for i = n-1:1
        x(i) = (y(i)-(gamma(i)*x(i+1)))/beta(i);
    end
end