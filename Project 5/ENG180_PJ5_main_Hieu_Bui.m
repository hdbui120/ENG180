clear, clc;
close all;

delta = [pi/10 pi/20 pi/40 pi/80];

[mesh,res1F] = finiteDiff(delta,pi/3,'forward');

figure('Name','Residual Plot')
loglog(mesh,res1F,'--')
xlabel('\pi/dx')
ylabel('Residual')
legend('Forward')






function [mesh,residual] = finiteDiff(deltax,x,style)
    
    % given function
    f = @(x) sin(x).^2;

    switch style
        case 'forward'
            mesh = pi./deltax;
    
            exactDiff = @(x) 2.*cos(x).*sin(x);
            dfdx = exactDiff(x);
            
            % First order forward
            diffF = (f(x+deltax)-f(x))./deltax;
            residual = abs(dfdx-diffF);
    end
end
