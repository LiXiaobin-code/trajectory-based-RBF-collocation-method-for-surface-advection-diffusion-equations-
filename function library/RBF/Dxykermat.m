function mat=Dxykermat(X, Y,coeff)
% by KCC
% creates D_{xy} of kernel matrix (3D only)
% for two point sets X and Y,Z
% The result is a matrix of size nx times ny
% coeff(1)Dxy+coeff(2)Dxz+coeff(3)Dyz
global RBFscale
[nx dx]=size(X);
[ny dy]=size(Y);
if dx~=dy
    error('Unequal space dimension for laplacekermat arguments');
end
% s=RBFscale.^2.*DistanceMatrixSquare(X,Y);


% 
% fmat1=frbf(DistanceMatrixSquare(X,Y)*diag(RBFscale.^2),1)*diag(RBFscale.^2);
fmat2=frbf(DistanceMatrixSquare(X,Y)*diag(RBFscale.^2),2)*diag(RBFscale.^4);
mat=zeros(nx,ny,dx);
if coeff(1)~=0
for dim=1:2
    mat(:,:,dim)=4*fmat2.*(repmat(X(:,dim),1,ny)-repmat(Y(:,dim)',nx,1)).*(repmat(X(:,3-dim),1,ny)-repmat(Y(:,3-dim)',nx,1));
end

elseif coeff(2)~=0

for dim=1:2:3
    mat(:,:,dim)=4*fmat2.*(repmat(X(:,dim),1,ny)-repmat(Y(:,dim)',nx,1)).*(repmat(X(:,end+1-dim),1,ny)-repmat(Y(:,end+1-dim)',nx,1));
end

else 
for dim=2:3
    mat(:,:,dim)=4*fmat2.*(repmat(X(:,dim),1,ny)-repmat(Y(:,dim)',nx,1)).*(repmat(X(:,end+2-dim),1,ny)-repmat(Y(:,end+2-dim)',nx,1));
end

end