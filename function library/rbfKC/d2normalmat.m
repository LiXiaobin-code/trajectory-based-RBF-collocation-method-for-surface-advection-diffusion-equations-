function mat=d2normalmat(X,N,Z)



[~,d]=size(X);



D2mat=D2allkermat(X,Z);

switch d
    
    case 2
   mat=diag(N(:,1).^2)*D2mat(:,:,1)+ diag(N(:,2).^2)*D2mat(:,:,2)+  2*diag(N(:,1).*N(:,2))*D2mat(:,:,3);
        
    case 3
   mat=diag(N(:,1).^2)*D2mat(:,:,1)+...
       diag(N(:,2).^2)*D2mat(:,:,2)+...
       diag(N(:,3).^2)*D2mat(:,:,3)+...
       2*diag(N(:,1).*N(:,2))*D2mat(:,:,4)+...
       2*diag(N(:,2).*N(:,3))*D2mat(:,:,5)+...
       2*diag(N(:,1).*N(:,3))*D2mat(:,:,6);     
        
        
end