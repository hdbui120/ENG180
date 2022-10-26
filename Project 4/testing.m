clc, clear; 
close all;


jacob = 'Jacobi';
gauss = 'GSeidel';
sor = 'SOR';
eqn = @(x,n) cos(pi.*x./(n+1));
[clk1,iterations1,res1,cor1] = iterationEfficiency(10,1.2,jacob,eqn);
[clk2,iterations2,res2,cor2] = iterationEfficiency(10,1.2,gauss,eqn);
[clk3,iterations3,res3,cor3] = iterationEfficiency(10,1.2,sor,eqn);
figure("Name",'Problem 1 ')
loglog(iterations1,res1)
hold on
loglog(iterations2,res2)
hold on
loglog(iterations3,res3)
hold on
title('Residuals vs Iterations')
legend1=sprintf('Jacobi      %fs ',clk1);
legend2=sprintf('Gauss S.    %fs',clk2);
legend3=sprintf('SOR         %fs',clk3);
xlabel('Iterations')
ylabel('Residuals')
legend(legend1,legend2,legend3,"Location",'southwest')


function [clk,iter,resValues, corValues] = iterationEfficiency(n,guess,method,eqn) 
    a = ones(n-1,1); 
    b = -2.*ones(n,1);
    c = ones(n-1,1);
    matrixA = diag(a,-1)+diag(b)+diag(c,1);
    iterations = 0;
    resid = 1;
    tolerance = 1e-15;
    f = eqn(1:n,n);
    xNew = guess.*ones(n,1);
    
    switch method
        case 'Jacobi'
            start = tic;
            while resid > tolerance
                xOld = xNew;
                for i=1:n
                    xNew(i) = xOld(i)+1/matrixA(i,i)*(f(i)-matrixA(i,:)*xOld(:));
                end
                resid = max(abs(matrixA*xNew-f'));
                cor = max(abs(xNew-xOld));
                iterations = iterations + 1;
                resValues(iterations) = resid;
                corValues(iterations) = cor;
            end
            clk = toc(start);
            iter = linspace(1,iterations,iterations);
        case 'GSeidel'
            start = tic;
            while resid > tolerance
                xOld = xNew;
                for i = 1:n
                    xNew(i) = xOld(i)+(f(i)-matrixA(i,1:i-1)*xNew(1:i-1)-matrixA(i,i:n)*xOld(i:n))/matrixA(i,i);
                end
                resid = max(abs(matrixA*xNew-f'));
                cor = max(abs(xNew-xOld));
                iterations = iterations + 1;
                resValues(iterations) = resid;
                corValues(iterations) = cor;
            end
            clk = toc(start);
            iter = linspace(1,iterations,iterations);
        case 'SOR'
            omega = 1.7;
            start = tic;
            while resid > tolerance
                xOld = xNew;
                for i = 1:n
                    xGS = xOld(i)+(f(i)-matrixA(i,1:i-1)*xNew(1:i-1)-matrixA(i,i:n)*xOld(i:n))/matrixA(i,i);
                    xNew(i) = omega*(xGS-xOld(i))+xOld(i);
                end
                resid = max(abs(matrixA*xNew-f'));
                cor = max(abs(xNew-xOld));
                iterations = iterations + 1;
                resValues(iterations) = resid;
                corValues(iterations) = cor;
            end
            clk = toc(start);
            iter = linspace(1,iterations,iterations);
    end
end