function D2=D2normalkermat(X,NX,Y)
% by KCC
% creates D_{xi xi} of kernel matrix
% for two point sets X and Y,
%2D: D(:,1)=x1x1;  D(:,2)=x2x2; D(:,2+1)=x1x2; D(:,2+2)=x2x1;
%3D: D(:,1)=x1x1;  D(:,2)=x2x2; D(:,3)=x3x3; D(:,3+1)=x1x2; D(:,3+2)=x2x3; D(:,3+3)=x3x1;
% The result is a matrix of size nx times ny
global RBFscale
[nx dx]=size(X);
[ny dy]=size(Y);
if dx~=dy
    error('Unequal space dimension for arguments');
end
% s=RBFscale.^2.*DistanceMatrixSquare(X,Y);



fmat1=frbf(DistanceMatrixSquare(X,Y)*diag(RBFscale.^2),1)*diag(RBFscale.^2);
fmat2=frbf(DistanceMatrixSquare(X,Y)*diag(RBFscale.^2),2)*diag(RBFscale.^4);
D2=zeros(nx,ny,dx+nchoosek(dx,2));
for dim=1:dx
    D2(:,:,dim)=4*fmat2.*(repmat(X(:,dim),1,ny)-repmat(Y(:,dim)',nx,1)).^2+2*fmat1;
end

for dim=dx+1:dx+nchoosek(dx,2)-1
    D2(:,:,dim)=4*fmat2.*(repmat(X(:,dim-dx),1,ny)-repmat(Y(:,dim-dx)',nx,1)).*(repmat(X(:,dim+1-dx),1,ny)-repmat(Y(:,dim+1-dx)',nx,1));
end

D2(:,:,dx+nchoosek(dx,2))=4*fmat2.*(repmat(X(:,dx),1,ny)-repmat(Y(:,dx)',nx,1)).*(repmat(X(:,1),1,ny)-repmat(Y(:,1)',nx,1));



switch dx
    
    case 2
        D2=diag(NX(:,1).^2)*D2(:,:,1)+diag(NX(:,2).^2)*D2(:,:,2)+2*diag(NX(:,1).*NX(:,2))*D2(:,:,3);    
        
    case 3
        D2=diag(NX(:,1).^2)*D2(:,:,1)+diag(NX(:,2).^2)*D2(:,:,2)+diag(NX(:,3).^2)*D2(:,:,3)+...
            2*(diag(NX(:,1).*NX(:,2))*D2(:,:,4)+diag(NX(:,2).*NX(:,3))*D2(:,:,5)+diag(NX(:,1).*NX(:,3))*D2(:,:,6)); 
        
end
      

