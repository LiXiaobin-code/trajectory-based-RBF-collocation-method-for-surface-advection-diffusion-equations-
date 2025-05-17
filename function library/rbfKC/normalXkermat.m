function mat=normalXkermat(X,NX,Y)
%Modified by KC CHEUNG, works on shape parameters epsilon
% create kernel matrices
% for two point sets X and Y
% corresponding to the normals NX wrt. the X varibale.
% The result is a matrix of size nx time ny.

global RBFscale
[nx dx]=size(X);
[ny dy]=size(Y);
if dx~=dy
    error('Unequal space dimension for normalXkermat arguments');
end
fmat=frbf(RBFscale^2.*DistanceMatrixSquare(X,Y),1)*RBFscale^2;

mat=2*fmat.*(repmat(diag(NX*X'),1,ny)-(NX*Y') );