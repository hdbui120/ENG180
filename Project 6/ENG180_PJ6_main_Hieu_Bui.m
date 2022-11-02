clc, close all;
clear;

n = 101;

%making the grid
x = linspace(-1,1,n);
delta = 2/n;
f = -4;
eps = .1;

iMinusOneCoeff = f*delta/2*eps + 1; 
iPlusOneCoeff = f*delta/2*eps - 1;
iCoeff = 2;

a = iMinusOneCoeff.*ones(n,1); b = iCoeff.*ones(n,1); c = iPlusOneCoeff.*ones(n,1);

d = zeros(n,1); d(1)=2; d(n) = 2;

%Solving tridiagonal system
u = THOMAS3(a,b,c,d,n);


% supporting function
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
        x(i) = (dbar(i)-(bbar(i)*x(i+1)))/cbar(i);
    end
end