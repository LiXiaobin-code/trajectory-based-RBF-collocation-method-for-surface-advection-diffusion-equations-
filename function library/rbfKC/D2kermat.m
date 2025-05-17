function mat=D2kermat(X, Y)
% by KCC
% creates D_{xi xi} of kernel matrix
% for two point sets X and Y,
% The result is a matrix of size nx times ny
global RBFscale
[nx dx]=size(X);
[ny dy]=size(Y);
if dx~=dy
    error('Unequal space dimension for laplacekermat arguments');
end
% s=RBFscale.^2.*DistanceMatrixSquare(X,Y);



fmat1=frbf(DistanceMatrixSquare(X,Y)*diag(RBFscale.^2),1)*diag(RBFscale.^2);
fmat2=frbf(DistanceMatrixSquare(X,Y)*diag(RBFscale.^2),2)*diag(RBFscale.^4);
mat=zeros(nx,ny,dx);
for dim=1:dx
    mat(:,:,dim)=4*fmat2.*(repmat(X(:,dim),1,ny)-repmat(Y(:,dim)',nx,1)).^2+2*fmat1;
end

