function normal = geomNormals(a,x,nref)
% This function calculates the normal and the curvature(s) of a polynomial
% at a given point
%
% Input:
%   x   : The point to calculate the geometric quantities.
%   a   : Polynomial coefficients
%   nref: Reference vector (for normal orientation)
%
% Output:
%   N   : The unit normal vector
%   curv: The curvature(s)
dim = length(x);

daX = ls_polyder(a,'x',dim);
    
if dim == 1
    normal = [ls_polyval(daX,x) -1];
    normal = normal./norm(normal,2);
    normal = rotate(normal,nref,2);
    normal = normal*sign(dot(normal,nref));
else
    daY = ls_polyder(a,'y',dim);
    normal = [ls_polyval(daX,x) ls_polyval(daY,x) -1];
    normal = normal./norm(normal,2);
    normal = rotate(normal,nref,2);
    normal = normal*sign(dot(normal,nref));
end