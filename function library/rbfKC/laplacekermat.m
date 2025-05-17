function mat=laplacekermat(X, Y)
% creates Laplacian of kernel matrix
% for two point sets X and Y,
% The result is a matrix of size nx times ny
global RBFscale
[nx dx]=size(X);
[ny dy]=size(Y);
if dx~=dy
    error('Unequal space dimension for laplacekermat arguments');
end
s=RBFscale^2*DistanceMatrixSquare(X,Y);
mat=(dx*frbf(s,1)+2*s.*frbf(s,2))*2*RBFscale^2;
