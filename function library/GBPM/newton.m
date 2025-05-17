function xmin = newton(a,xmin)
% This function implements a Newton-Raphson root finder algorithm to solve
% the minimization problem
% 
% d(x) = 1/2*||(x,f(x))||^2
%
% Input:
%   xmin : The initial guess.
%   a    : Polynomial coefficients
%
% Output:
%   xmin : The minimizer of the polynomial.

n = length(xmin);
dim = length(xmin);

daX = ls_polyder(a,'x',dim);
d2aXX = ls_polyder(daX,'x',dim);
if n == 1
    counter=0; err=1;
    while err>1e-8 && counter<100
        xminPrev = xmin;
        F = xmin + ls_polyval(a,xmin)*ls_polyval(daX,xmin);
        J = 1 + ls_polyval(daX,xmin).^2 + ls_polyval(a,xmin)*ls_polyval(d2aXX,xmin);
        xmin = xmin - J\F;
        err = abs(xmin - xminPrev)/max([abs(xmin),abs(xminPrev),1e-10]);
        counter = counter + 1;
    end
elseif n == 2
    daY = ls_polyder(a,'y',dim);
    d2aYY = ls_polyder(daY,'y',dim);
    d2aXY = ls_polyder(daX,'y',dim);
    
    counter=0; err1=1; err2=1;
    while err1>1e-8 && err2>1e-8 && counter<100
        xminPrev = xmin;
        F = xmin + [ls_polyval(a,xmin')*ls_polyval(daX,xmin');ls_polyval(a,xmin')*ls_polyval(daY,xmin')];
        J = [1 + ls_polyval(daX,xmin').^2 + ls_polyval(a,xmin')*ls_polyval(d2aXX,xmin'),...
            ls_polyval(daY,xmin')*ls_polyval(daX,xmin') + ls_polyval(a,xmin')*ls_polyval(d2aXY,xmin');...
            ls_polyval(daY,xmin')*ls_polyval(daX,xmin') + ls_polyval(a,xmin')*ls_polyval(d2aXY,xmin'),...
            1 + ls_polyval(daY,xmin').^2 + ls_polyval(a,xmin')*ls_polyval(d2aYY,xmin')];
        xmin = xmin - J\F;
        err1 = abs(xmin(1) - xminPrev(1))/max([abs(xmin(1)),abs(xminPrev(1)),1e-10]);
        err2 = abs(xmin(2) - xminPrev(2))/max([abs(xmin(2)),abs(xminPrev(2)),1e-10]);
        counter = counter + 1;
    end    
    xmin = xmin';
else
    error('Newton''s method not implemented for higher dimensions') 
end