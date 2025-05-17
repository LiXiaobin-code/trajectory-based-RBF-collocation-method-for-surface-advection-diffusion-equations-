function xout = rotate(xin,nref,par)
% This function rotates the given points in the local coordinates given by
% the reference normal vector using the Householder reflection method.
%
% Input:
%   xin : The point(s) to be rotated.
%   nref: The reference normal vector for the rotation.
%   par : 1 to rotate to local coordinates, 2 to rotate to Cartesian
%       coordinates.
%
% Output:
%   xout: The rotated point(s).

dim = size(nref,2);
nref = nref/norm(nref);

if dim == 2
    A = [nref(2) nref(1);-nref(1) nref(2)]';
elseif dim == 3
    % Householder matrix method to get the rotation matrix
    I = eye(3);
    u = (nref-I(end,:))/(2*sqrt((1/2)*(1+abs(nref(3)))));
    if norm(u)~=1 && norm(u)>1e-15
        u = u/norm(u);
    end
    % Rotation matrix
    A = I-2*(u'*u);
else
    error('Not implemented for higher dimensions')
end

if par == 1 % Rotate to local coordinates
    xout = (A*xin')';
else % Rotate to Cartesian coordinates
    xout = (A'*xin')';
end