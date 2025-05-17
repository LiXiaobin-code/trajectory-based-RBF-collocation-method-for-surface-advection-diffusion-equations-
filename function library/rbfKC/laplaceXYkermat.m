function mat=laplaceXYkermat(X, Y)
% creates Laplacian of kernel matrix
% for two point sets X and Y,
% the Laplacian applied to BOTH arguments.
% The result is a matrix of size nx times ny
global RBFscale
[nx dx]=size(X);
[ny dy]=size(Y);
if dx~=dy
    error('Unequal space dimension for laplaceXYkermat arguments');
end
s=RBFscale^2*DistanceMatrixSquare(X,Y);
mat=(dx*(dx+2)*frbf(s,2)+4*(dx+2)*s.*frbf(s,3)+4*s.^2.*frbf(s,4))*4*RBFscale^4;
