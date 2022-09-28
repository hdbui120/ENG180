clear, clc; 
close all;

n = 4;
a = 2.*ones(1,n); % a(1) = 1
b = -5.*ones(1,n);
c = 2.*ones(1,n);
d = [-3,-1,-1,-3];
%zBefore = diag(a,-1) + diag(b) + diag(c,1);
%zBefore = [zBefore, d']
z = THOMAS3(a,b,c,d,n);

function x = THOMAS3(a,b,c,d,n)

    %initialize first row
    b(1) = b(1);
    c(1) = c(1);
    
    %making upper triangle
    for i = 2:n
        multiplier = a(i)./b(i-1);
        a(i) = a(i) - b(i-1).*multiplier;
        b(i) = b(i) - c(i-1).*multiplier;
        c(i) = c(i);
        d(i) = d(i) - d(i-1).*multiplier;
    end
    
    %initialize x of size n
    x = ones(1,n);

    %initialize end condition
    x(1) = d(n)/b(n);

    % Upward substitution AKA zip it up
    for i = n-1:1
        x(i) = (d(i)-(b(i)*x(i+1)))/c(i);
    end

end