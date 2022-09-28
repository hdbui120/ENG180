clear, clc; 
close all;

%problem 3

%3A
%declare and initialize variables
a
b
c
d
e
f
n

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
        bbar(i) = b(i)-multiplier*dbar(i-2);
        cbar(i) = c(i)-multiplier*ebar(i-2);
        dbar(i) = d(i);
        ebar(i) = e(i);
        fbar(i) = f(i)-multiplier*fbar(i-2);
    end

%     downward elimination for b 
    for i = 3:n
        multiplier2 = bbar(i)/cbar(i-1);
        cbar(i) = c(i)-multiplier2*dbar(i-1);
        dbar(i) = dbar(i)-multiplier2*ebar(i-1);
        ebar(i) = e(i);
        fbar(i) = fbar(i)-multiplier2*fbar(i-1);
    end    

%     preallocate x

%     initial conditions for upward substitution
    
end