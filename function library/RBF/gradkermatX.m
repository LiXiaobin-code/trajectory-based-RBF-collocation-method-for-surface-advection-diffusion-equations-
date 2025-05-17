function mat=gradkermatX(X,Y,scale)
[nx dx]=size(X);
[ny dy]=size(Y);
if dx~=dy
    error('Unequal space dimension for gradkermatX arguments');
end
fmat=frbf(DistanceMatrixSquare(X,Y)*diag(scale.^2),1)*diag(scale.^2);
mat=zeros(nx,ny,dx);
for dim=1:dx
    mat(:,:,dim)=2*fmat.*(repmat(X(:,dim),1,ny)-repmat(Y(:,dim)',nx,1));
end