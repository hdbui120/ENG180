clear, clc; 
close all;

%declaring and initializing variables
n = 4;
a = 2.*ones(1,n); % a(1) = 1
b = -5.*ones(1,n);
c = 2.*ones(1,n);
d = [-3,-1,-1,-3];
z = THOMAS3(a,b,c,d,n);

%problem 1
function x = THOMAS3(a,b,c,d,n)

    %initialize first row
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
    for i = n-1:1
        x(i) = (dbar(i)-(bbar(i)*x(i+1)))/cbar(i);
    end
end
% end problem 1

% begin problem 2
% 2x2 inverse function
% function to find the inverse of 2x2 matrix
function in = invert(m)
    determ = (m(1,1).*m(2,2))-(m(1,2).*m(2,1));
    rearrange = [m(2,2), -1.*m(1,2); -1.*m(2,1), m(1,1)];
    in = (1/determ).*rearrange;
end

function y = btrid(a,b,c,d,n)
    
    %initialize first row
    bbar{1} = b{1};
    cbar{1} = c{1};
    dbar{1} = d{1};

    %downward elimination
    for i = 2:n
        multiplier = a{i}*invert(bbar{i-1});
        abar{i} = a{i} - multiplier*bbar{i-1};
        bbar{i} = b{i} - multiplier*cbar{i-1};
        cbar{i} = c{i};
        dbar{i} = d{i} - multiplier*dbar{i-1};
    end
    
    %initialize cell array with matrix of size 2x2
    y = cell(n,1);
    for i = 1:n
        y{i} = [1;1];
    end
    
    %initilize end condition for upward substitution
    y{n} = invert(bbar{n})*dbar{n};
    
    %upward substitution
    for i = n-1:1
        y{i} = invert(cbar{i})*(dbar{i}-(bbar{i}*y{i+1}));
    end
end
% end problem2
