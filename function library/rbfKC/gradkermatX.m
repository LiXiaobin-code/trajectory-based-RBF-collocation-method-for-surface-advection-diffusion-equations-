% function mat=gradkermatX(X,Y)
% % create kernel matrices
% % for two points sets X and Y
% % corresponding to full gradient wrt. the X variable.
% % The result is a 3-dimensional matrix of size nx times nt times dx=dy,
% 
% global RBFscale
% [nx dx]=size(X);
% [ny dy]=size(Y);
% if dx~=dy
%     error('Unequal space dimension for gradkermatX arguments');
% end
% %fmat=frbf(RBFscale.^2.*DistanceMatrixSquare(X,Y),1).*RBFscale.^2;
% fmat=frbf(DistanceMatrixSquare(X,Y)*diag(RBFscale.^2),1)*diag(RBFscale.^2);
% mat=zeros(nx,ny,dx);
% for dim=1:dx
%     mat(:,:,dim)=2*fmat.*(repmat(X(:,dim),1,ny)-repmat(Y(:,dim)',nx,1));
% end

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