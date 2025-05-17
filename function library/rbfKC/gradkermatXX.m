function mat=gradkermatXX(X,Y)
% create kernel matrices
% for two points sets X and Y
% corresponding to full gradient wrt. the X variable.
% The result is a 3-dimensional matrix of size nx times nt times dx=dy,

global RBFscale
[nx dx]=size(X);
[ny dy]=size(Y);
if dx~=dy
    error('Unequal space dimension for gradkermatX arguments');
end
fmat1=frbf(RBFscale^2*DistanceMatrixSquare(X,Y),1)*RBFscale^2;
fmat2=frbf(RBFscale^2*DistanceMatrixSquare(X,Y),2)*RBFscale^4;
mat=zeros(nx,ny,dx);
for dim=1:dx
    mat(:,:,dim)=2*fmat1+4*fmat2.*(repmat(X(:,dim),1,ny)-repmat(Y(:,dim)',nx,1)).^2;
end