clear, clc;
close all;

delta = [pi/10 pi/20 pi/40 pi/80];

% Compute 1st order derivative residuals
res1F1 = finiteDiff(delta,pi/3,'forward1');
res1B1 = finiteDiff(delta,pi/3,'backward1');
res1C2 = finiteDiff(delta,pi/3,'central1');
res1F2 = finiteDiff(delta,pi/3,'forward2');
res1B2 = finiteDiff(delta,pi/3,'backward2');
resCubic = finiteDiff(delta,pi/3,'cubic');

%Compute 2nd order derivative residuals
res2F1 = finiteDiff2(delta,pi/3,'forward1');
res2B1 = finiteDiff2(delta,pi/3,'backward1');
res2C2 = finiteDiff2(delta,pi/3,'central');
res2F2 = finiteDiff2(delta,pi/3,'forward2');
res2B2 = finiteDiff2(delta,pi/3,'backward2');
resCubic2 = finiteDiff2(delta,pi/3,'cubic');

% Compute 3rd order derivative residuals
res3F1 = finiteDiff3(delta,pi/3,'forward1');
res3B1 = finiteDiff3(delta,pi/3,'backward1');
res3C2 = finiteDiff3(delta,pi/3,'central');
res3F2 = finiteDiff3(delta,pi/3,'forward2');
res3B2 = finiteDiff3(delta,pi/3,'backward2');

% Residual for fourth order derivative
res4C2 = finiteDiff4(delta,pi/3,'central');
res4F1 = finiteDiff4(delta,pi/3,'forward');
res4B1 = finiteDiff4(delta,pi/3,'backward');

figure('Name','Residual Plot 1st Order')
loglog(res1F1,'Marker','o')
hold on
loglog(res1B1,'Marker','*')
hold on
loglog(res1C2,'Marker','diamond')
hold on
loglog(res1F2,'Marker','<')
hold on
loglog(res1B2,'Marker','+')
hold on
loglog(resCubic,'Marker','^')
title('Problem 1a: First Derivative of f(x) = $sin^2(x)$ at $x=\frac{\pi}{3}$','Interpreter','latex')
xlabel('\pi/dx')
ylabel('Residuals')
legend('Forward','Backward','Central','2nd Order Forward','2nd Order Backward','Cubic Spline')

figure('Name','Residual Plot 2nd Order')
loglog(res2F1,'Marker','o')
hold on
loglog(res2B1,'Marker','*')
hold on
loglog(res2C2,'Marker','diamond')
hold on
loglog(res2F2,'Marker','<')
hold on
loglog(res2B2,'Marker','+')
hold on
loglog(resCubic2,'Marker','^')
title('Problem 1a: Second Derivative of f(x) = $sin^2(x)$ at $x=\frac{\pi}{3}$','Interpreter','latex')
xlabel('\pi/dx')
ylabel('Residuals')
legend('Forward','Backward','Central','2nd Order Forward','2nd Order Backward','Cubic Spline')

figure('Name','Residual Plot 3rd Order')
loglog(res3F1,'Marker','o')
hold on
loglog(res3B1,'Marker','*')
hold on
loglog(res3C2,'Marker','diamond')
hold on
loglog(res3F2,'Marker','<')
hold on
loglog(res3B2,'Marker','+')
title('Problem 1a: Third Derivative of f(x) = $sin^2(x)$ at $x=\frac{\pi}{3}$','Interpreter','latex')
xlabel('\pi/dx')
ylabel('Residuals')
legend('Forward','Backward','Central','2nd Order Forward','2nd Order Backward')

figure('Name','Residual Plot 4th Order')
loglog(res4C2,'Marker','diamond')
hold on
loglog(res4F1,'Marker','o')
hold on
loglog(res4B1,'Marker','*')
title('Problem 1a: Fourth Derivative of f(x) = $sin^2(x)$ at $x=\frac{\pi}{3}$','Interpreter','latex')
xlabel('\pi/dx')
ylabel('Residuals')
legend('Central','Forward','Backward')

% main function to find 1st order derivative for problem 1A
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
        case 'cubic'
            % First derivative cubic spline
            k1=exactDiff(x-deltax);
            k3 = exactDiff(x+deltax);
            diffF = 3.*(f(x+deltax)-f(x-deltax))./(4.*deltax)-((k1+k3)./4);
            residual = abs(dfdx-diffF);
    end
end

% main function to find 2nd order derivative for problem 1A
function residual = finiteDiff2(delta,x,style)
    % given function
    f = @(x) sin(x).^2;
    exactDiff = @(x) 2.*(cos(x).^2-sin(x).^2);
    div1 = @(x) 2.*cos(x).*sin(x);
    dfdx = exactDiff(x);
    k1=div1(x-delta);
    k3 = div1(x+delta);
    ki = 3.*(f(x+delta)-f(x-delta))./(4.*delta)-((k1+k3)./4);

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
        case 'cubic'
            % First derivative cubic spline
            diffF = 3.*(f(x+delta)-f(x))./(delta.^2)-(k3+2.*ki)./delta;
            residual = abs(dfdx-diffF);
    end
end

% main function to find 3rd order derivative for problem 1A
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
            diffF = (-.5.*f(x-2.*delta)+f(x-delta)-f(x+delta)+.5.*f(x+2.*delta))./delta.^3;
            residual = abs(dfdx-diffF);
        case 'forward2'
            % Second order of accurate 
            diffF = ((-5/2).*f(x)+9.*f(x+delta)-12.*f(x+2.*delta)+7.*f(x+3.*delta)-(3/2).*f(x+4.*delta))./delta.^3;
            residual = abs(dfdx-diffF);
        case 'backward2'
            % Second order of accurate 
            diffF = (5/2.*f(x)-9.*f(x-delta)+12.*f(x-2.*delta)-7.*f(x-3.*delta)+3/2.*f(x-4.*delta))./delta.^3;
            residual = abs(dfdx-diffF);
    end
end

% main function to find 4th order derivative for problem 1A
function residual = finiteDiff4(delta,x,style)
    % given function
    f = @(x) sin(x).^2;
    exactDiff = @(x) -8.*(cos(x).^2-sin(x).^2);
    dfdx = exactDiff(x);
    
    switch style
        case 'central'
            % Second order of accurate central
            diffF = (f(x-2.*delta)-4.*f(x-delta)+6.*f(x)-4.*f(x+delta)+f(x+2.*delta))./delta.^4;
            residual = abs(dfdx-diffF);
        case 'forward'
            % First order of accurate
            diffF = (f(x)-4.*f(x+delta)+6.*f(x+2.*delta)-4.*f(x+3.*delta)+f(x+4.*delta))./delta.^4;
            residual = abs(dfdx-diffF);
        case 'backward'
            % Second order of accurate central
            diffF = (f(x)-4.*f(x-delta)+6.*f(x-2.*delta)-4.*f(x-3.*delta)+f(x-4.*delta))./delta.^4;
            residual = abs(dfdx-diffF);
    end
end





