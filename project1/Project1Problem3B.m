clear, clc;
close all;

%3B
n = 4;
for i = 1:n
    a{i} = [1 -2; 0 1];
    b{i} = [-5 -1; -2 -5];
    c{i} = [1 0; -1 1];
end
f = {[-5;-7];[-6;-6];[-6;-6];[-7;-6]};

z = btrid(a,b,c,f,n);
celldisp(z)

% 2x2 inverse function
% function to find the inverse of 2x2 matrix
function in = invert(m)
    determ = (m(1,1).*m(2,2))-(m(1,2).*m(2,1));
    rearrange = [m(2,2), -1.*m(1,2); -1.*m(2,1), m(1,1)];
    in = (1/determ).*rearrange;
end

function y = btrid(a,b,c,f,n)
    
    %initialize first row
    bbar{1} = b{1};
    cbar{1} = c{1};
    dbar{1} = f{1};

    %downward elimination
    for i = 2:n
        multiplier = a{i}*invert(bbar{i-1});
        abar{i} = a{i} - multiplier*bbar{i-1};
        bbar{i} = b{i} - multiplier*cbar{i-1};
        cbar{i} = c{i};
        dbar{i} = f{i} - multiplier*dbar{i-1};
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