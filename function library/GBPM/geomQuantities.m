function [normal,curv,tang] = geomQuantities(a,x,nref)
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
d2aXX = ls_polyder(daX,'x',dim);
    
if dim == 1
    normal = [ls_polyval(daX,x) -1];
    normal = normal./norm(normal,2);
    normal = rotate(normal,nref,2);
    normal = normal*sign(dot(normal,nref));
    
    tang = [-normal(2) normal(1)];

    curv = ls_polyval(d2aXX,x)/sqrt(1+ls_polyval(daX,x)^2)^3;
else
    daY = ls_polyder(a,'y',dim);
    d2aYY = ls_polyder(daY,'y',dim);
    d2aXY = ls_polyder(daX,'y',dim);
    
    normal = [ls_polyval(daX,x) ls_polyval(daY,x) -1];
    normal = normal./norm(normal,2);
    normal = rotate(normal,nref,2);
    normal = normal*sign(dot(normal,nref));
    
    % Householder matrix method to get the rotation matrix
    u = [normal(1:2) max(normal(3)+norm(normal),normal(3)-norm(normal))];
    I = eye(3);
%     u = (normal-I(end,:))/(2*sqrt((1/2)*(1+abs(normal(3)))));
%     if norm(u)~=1 && norm(u)>1e-15
%         u = u/norm(u);
%     end
    % Tangential vector matrix (tangent vectors as rows)
    tang = I-2*(u'*u)/(norm(u)^2);
    
    fx = ls_polyval(daX,x); fy = ls_polyval(daY,x);
    fxx = ls_polyval(d2aXX,x); fxy = ls_polyval(d2aXY,x);
    fyx = ls_polyval(d2aXY,x); fyy = ls_polyval(d2aYY,x);
    K = (fxx*fyy-fxy*fyx)/((1+fx^2+fy^2)^2);
    H = ((1+fx^2)*fyy-fx*fy*fxy-fy*fx*fyx+(1+fy^2)*fxx)/((1+fx^2+fy^2)^(3/2));
    curv = [H K];    
end