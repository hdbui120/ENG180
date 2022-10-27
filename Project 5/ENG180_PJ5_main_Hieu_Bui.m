clear, clc;
close all;

delta = [pi/10 pi/20 pi/40 pi/80];

% Compute 1st order derivative residuals
res1F1 = finiteDiff(delta,pi/3,'forward1');
res1B1 = finiteDiff(delta,pi/3,'backward1');
res1C1 = finiteDiff(delta,pi/3,'central1');
res1F2 = finiteDiff(delta,pi/3,'forward2');
res1B2 = finiteDiff(delta,pi/3,'backward2');

%Compute 2nd order derivative residuals
res2F1 = finiteDiff2(delta,pi/3,'forward1');
res2B1 = finiteDiff2(delta,pi/3,'backward1');
res2C1 = finiteDiff2(delta,pi/3,'central');
res2F2 = finiteDiff2(delta,pi/3,'forward2');
res2B2 = finiteDiff2(delta,pi/3,'backward2');

% Compute 3rd order derivative residuals
res3F1 = finiteDiff3(delta,pi/3,'forward1');
res3B1 = finiteDiff3(delta,pi/3,'backward1');
res3C1 = finiteDiff3(delta,pi/3,'central');
res3F2 = finiteDiff3(delta,pi/3,'forward2');
res3B2 = finiteDiff3(delta,pi/3,'backward2');

figure('Name','Residual Plot 1st Order')
loglog(res1F1,'Marker','o')
hold on
loglog(res1B1,'Marker','*')
hold on
loglog(res1C1,'Marker','diamond')
hold on
loglog(res1F2,'Marker','<')
hold on
loglog(res1B2,'Marker','+')
title('Problem 1a: First Derivative of f(x) = $sin^2(x)$ at $x=\frac{\pi}{3}$','Interpreter','latex')
xlabel('\pi/dx')
ylabel('Residuals')
legend('Forward','Backward','Central','2nd Order Forward','2nd Order Backward')

figure('Name','Residual Plot 2nd Order')
loglog(res2F1,'Marker','o')
hold on
loglog(res2B1,'Marker','*')
hold on
loglog(res2C1,'Marker','diamond')
hold on
loglog(res2F2,'Marker','<')
hold on
loglog(res2B2,'Marker','+')
title('Problem 1a: Second Derivative of f(x) = $sin^2(x)$ at $x=\frac{\pi}{3}$','Interpreter','latex')
xlabel('\pi/dx')
ylabel('Residuals')
legend('Forward','Backward','Central','2nd Order Forward','2nd Order Backward')

figure('Name','Residual Plot 3rd Order')
loglog(res3F1,'Marker','o')
hold on
loglog(res3B1,'Marker','*')
hold on
loglog(res3C1,'Marker','diamond')
hold on
loglog(res3F2,'Marker','<')
hold on
loglog(res3B2,'Marker','+')
title('Problem 1a: Third Derivative of f(x) = $sin^2(x)$ at $x=\frac{\pi}{3}$','Interpreter','latex')
xlabel('\pi/dx')
ylabel('Residuals')
legend('Forward','Backward','Central','2nd Order Forward','2nd Order Backward')

function residual = finiteDiff(deltax,x,style)    
    % given function
    f = @(x) sin(x).^2;
    exactDiff = @(x) 2.*cos(x).*sin(x);
    dfdx = exactDiff(x);
    switch style
        case 'forward1'       
            % First order of accurate forward
            diffF = (f(x+deltax)-f(x))./deltax;
            residual = abs(dfdx-diffF);

        case 'backward1'
            % First order of accurate backward
            diffF = (f(x)-f(x-deltax))./deltax;
            residual = abs(dfdx-diffF);
            
        case 'central1'
            % Second order of accurate central
            diffF = (f(x+deltax)-f(x-deltax))./(2.*deltax);
            residual = abs(dfdx-diffF);
        case 'forward2'
            % Second order of accurate forward
            diffF = (-3.*f(x)+4.*f(x+deltax)-f(x+(2.*deltax)))./(2.*deltax);
            residual = abs(dfdx-diffF);
        case 'backward2'
            % Second order of accurate backward
            diffF = (3.*f(x)-4.*f(x-deltax)+f(x-(2.*deltax)))./(2.*deltax);
            residual = abs(dfdx-diffF);
    end
end

function residual = finiteDiff2(delta,x,style)
    % given function
    f = @(x) sin(x).^2;
    exactDiff = @(x) 2.*(cos(x).^2-sin(x).^2);
    dfdx = exactDiff(x);
    
    switch style
        case 'forward1'
            % First order of accurate
            diffF = (f(x)-2.*f(x+delta)+f(x+2.*delta))./delta.^2;
            residual = abs(dfdx-diffF);
        case 'backward1'
            % First order of accurate
            diffF = (f(x)-2.*f(x-delta)+f(x-2.*delta))./delta.^2;
            residual = abs(dfdx-diffF);
        case 'central'
            % Second order of accurate central
            diffF = (f(x-delta)-2.*f(x)+f(x+delta))./delta.^2;
            residual = abs(dfdx-diffF);
        case 'forward2'
            % Second order of accurate 
            diffF = (2.*f(x)-5.*f(x+delta)+4.*f(x+2.*delta)-f(x+3.*delta))./delta.^2;
            residual = abs(dfdx-diffF);
        case 'backward2'
            % Second order of accurate 
            diffF = (2.*f(x)-5.*f(x-delta)+4.*f(x-2.*delta)-f(x-3.*delta))./delta.^2;
            residual = abs(dfdx-diffF);
    end
end

function residual = finiteDiff3(delta,x,style)
    % given function
    f = @(x) sin(x).^2;
    exactDiff = @(x) -8.*cos(x).*sin(x);
    dfdx = exactDiff(x);
    
    switch style
        case 'forward1'
            % First order of accurate
            diffF = (-f(x)+3.*f(x+delta)-3.*f(x+2.*delta)+f(x+3.*delta))./delta.^3;
            residual = abs(dfdx-diffF);
        case 'backward1'
            % First order of accurate
            diffF = (f(x)-3.*f(x-delta)+3.*f(x-2.*delta)-f(x-3.*delta))./delta.^3;
            residual = abs(dfdx-diffF);
        case 'central'
            % Second order of accurate central
            diffF = (-.5.*f(x-2.*delta)+f(x-delta)+f(x)-f(x+delta)+.5.*f(x+2.*delta))./delta.^2;
            residual = abs(dfdx-diffF);
        case 'forward2'
            % Second order of accurate 
            diffF = (-5/2.*f(x)+9.*f(x+delta)-12.*f(x+2.*delta)+7.*f(x+3.*delta)-3/2.*f(x+4.*delta))./delta.^2;
            residual = abs(dfdx-diffF);
        case 'backward2'
            % Second order of accurate 
            diffF = (5/2.*f(x)-9.*f(x-delta)+12.*f(x-2.*delta)-7.*f(x-3.*delta)+3/2.*f(x-4.*delta))./delta.^2;
            residual = abs(dfdx-diffF);
    end
end





