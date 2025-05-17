function a = ls_polyfit(x,y,deg)
% This function generates the least-squares polynomial fit in 1 and 2
% dimensions
%
% Input:
%   x    : The independent variables (Nx1 or Nx2)
%   y    : The dependent variable (Nx1)
%   deg  : The degree of the fitting polynomial
%
% Output:
%   a    : The coefficients of the least-squares polynomial
%
% The polynomials are of the form:
%   1D: y(x) = c + x + x^2 + ...
%   2D: z(x,y) = c + y + y^2 + ... + x + xy + xy^2 + ... 
%       (e.g. for 2nd order, z(x,y) = c + y + y^2 + x + xy + x^2)

% Ensuring corret input and least squares approach
[m,n] = size(x);
if n>m
    x = x';
    [~,n] = size(x);
end

if n == 1
    % 1D case is easy
    M = x.^(0:deg);
elseif n == 2
    % Identify the polynomial coefficients in 2D
    n_x = repelem(0:deg,(deg+1):-1:1);
    n_y = repmat(0:deg,1,deg+1) - repelem(0:deg,deg+1);
    n_y = n_y(n_y>=0);
    
    % Construct polynomial
    M = x(:,1).^(n_x).*x(:,2).^(n_y);
end

% Least-squares fit using the QR factorization
[Q,R,P] = qr(M,0);
a(P,:) = R\(Q'*y);