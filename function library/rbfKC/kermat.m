function mat=kermat(X, Y)
% creates kernel matrix
% for two point sets X and Y  
global RBFscale
[nx dx]=size(X);
[ny dy]=size(Y);
if dx~=dy
    error('Unequal space dimension for kermat arguments');
end
mat=frbf(RBFscale^2*DistanceMatrixSquare(X,Y),0);
