function f = ls_polyval(a,x)
% This function evaluates a polynomial that was constructed from 
% the ls_polyfit function at a point (or points) x
%
% Input:
%   a  : Polynomial coefficients (Nx1)
%   x  : Points to be evaluated (Mx1 or Mx2)
%
% Output:
%   f  : The evaluated function

n = length(x(1,:)); % Number of points to evaluate
N = length(a); % Number of coefficients

if n == 1
    % Reverse engineer degree 
    deg = N-1;
    % Construct polynomial matrix
    M = x.^(0:deg);
elseif n == 2
    % Reverse engineer degree 
    deg = round((-3+sqrt(1 + 8*N))/2);
    % Construct polynomial matrix
    n_x = repelem(0:deg,(deg+1):-1:1);
    n_y = repmat(0:deg,1,deg+1) - repelem(0:deg,deg+1);
    n_y = n_y(n_y>=0);
    M = x(:,1).^(n_x).*x(:,2).^(n_y);
else
    error('Not implemented for higher than 3 dimensions');
end

% Evaluate function
f = M*a;