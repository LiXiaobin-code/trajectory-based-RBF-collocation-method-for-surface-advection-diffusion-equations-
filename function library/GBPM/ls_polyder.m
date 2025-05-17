function aDer = ls_polyder(a,der,dim)
% This function gives the coefficients of the derivatives of a polynomial 
% that was constructed from the ls_polyfit function. This function can only
% calculate one derivative at a time
%
% Input:
%   a   : Polynomial coefficients (Nx1)
%   der : Derivative to be calculated ('x' or 'y')
%   dim : Dimensions of the polynomial
%
% Output:
%   f  : The evaluated function

% Number of polynomial coefficients
n = length(a);

if dim == 1
    % Easy to calculate the derivative in the 1D case
    aDer = a(2:n).*(1:(n-1))';
elseif dim == 2
    % Reverse engineer polynomial degree
    deg = round((-3+sqrt(1 + 8*n))/2);
    
    if der == 'x' % First derivative in x
        % Adjust polynomial degree
        v = repelem(0:deg,(deg+1):-1:1)';
        aDer = a.*v;
        % Find coefficients for the derivative
        aDer = aDer(v>0); 
    elseif der == 'y'
        % Adjust polynomial degree
        v = repmat(0:deg,1,deg+1)' - repelem(0:deg,deg+1)';
        v = v(v>=0);
        aDer = a.*v;
        % Find coefficients for the derivative
        aDer = aDer(v>0);                
    else
        error('Invalid option. Please select among the options for the derivative: "x", "y"');
    end
end