clear, clc; 
close all;

%problem 3

%3A
%declare and initialize variables
n = 8;
a = ones(1:n);
b = -2.*ones(1:n);
c = -5.*ones(1:n);
d = -1.*ones(1:n);
e = ones(1:n);
f = [-5 -7 -6 -6 -6 -6 -7 -6];

z = pent(a,b,c,d,e,f,n)

function x = pent(a,b,c,d,e,f,n)

%     initial condition 1
    cbar(1) = c(1);
    dbar(1) = d(1);
    ebar(1) = e(1);
    fbar(1) = f(1);

%     initial condition 2
    cbar(2) = c(2)-(b(2)/cbar(1)*dbar(1));
    dbar(2) = d(2)-(b(2)/cbar(1)*ebar(1));
    ebar(2) = e(2);
    fbar(2) = f(2)-(b(2)/cbar(1)*fbar(1));

%     downward elimination for a
    for i = 3:n
        multiplier = a(i)/cbar(i-2);
        abar(i) = a(i)-multiplier*cbar(i-2);
        bbar(i) = b(i)-multiplier*dbar(i-2);
        cbar(i) = c(i)-multiplier*ebar(i-2);
        dbar(i) = d(i);
        ebar(i) = e(i);
        fbar(i) = f(i)-multiplier*fbar(i-2);
    end
    
%     downward elimination for b 
    for i = 3:n
        multiplier2 = bbar(i)/cbar(i-1);
        bbar2(i) = bbar(i)-multiplier2*cbar(i-1);
        cbar2(i) = cbar(i)-multiplier2*dbar(i-1);
        dbar2(i) = dbar(i)-multiplier2*ebar(i-1);
        ebar2(i) = e(i);
        fbar2(i) = fbar(i)-multiplier2*fbar(i-1);
    end   

%     preallocate x
    x = ones(1,n);

%     initial conditions for upward substitution
    x(n) = fbar2(n)/cbar2(n);
    x(n-1) = (fbar2(n-1)-(dbar2(n-1)*x(n)))/cbar2(n-1);

%     upward substitution
    for i = n-2:1
        x(i) = (fbar2(i)-ebar2(i)*x(i+2)-dbar2(i)*x(i+1))/cbar2(i);
    end
end