function P=GetPowFun(X,Y,lambda,op)
% X: center point
% Y: local qneighbor point

global dim RBFscale RBFtype ;

%RBFscale=2;
lambda=lambda';

switch op
    
    
    
    case 'x'
    s=RBFscale^2*DistanceMatrixSquare(X,Y)';  
    sx=(repmat(X(:,1),size(Y,1),1)-Y(:,1));
    mat=2*RBFscale^2*sx.*frbf(s,1);  
    

        switch RBFtype
            case 'g'
    P=-2*RBFscale^2*frbf(0,1)-lambda'*mat;
            case 'ms'
    P=-2*RBFscale^2*frbf(0,1)-2*lambda'*mat+...
        lambda'*frbf(RBFscale^2*DistanceMatrixSquare(Y,Y),0)*lambda;                
                
                
        end
    
    
    case 'y'
        
    s=RBFscale^2*DistanceMatrixSquare(X,Y)';  
    sx=(repmat(X(:,2),size(Y,1),1)-Y(:,2));
    mat=2*RBFscale^2*sx.*frbf(s,1);    
   
    
        switch RBFtype
            case 'g'
    P=-2*RBFscale^2*frbf(0,1)-lambda'*mat;
            case 'ms'
    P=-2*RBFscale^2*frbf(0,1)-2*lambda'*mat+...
        lambda'*frbf(RBFscale^2*DistanceMatrixSquare(Y,Y),0)*lambda;                
                
                
        end
            
        
        
    case 'xx'
        
    s=RBFscale^2*DistanceMatrixSquare(X,Y)';
    
    sxx=(repmat(X(:,1),size(Y,1),1)-Y(:,1)).^2;
    
    mat=2*RBFscale^2*frbf(s,1)+4*RBFscale^4*sxx.*frbf(s,2);   
    
        switch RBFtype
            case 'g'
   P=12*RBFscale^4*frbf(0,2)-lambda'*mat;
            case 'ms'
    P=12*RBFscale^4*frbf(0,2)-2*lambda'*mat+...
        lambda'*frbf(RBFscale^2*DistanceMatrixSquare(Y,Y),0)*lambda;                
                
                
        end  
    
    
    case 'yy'
    s=RBFscale^2*DistanceMatrixSquare(X,Y)';
    
    syy=(repmat(X(:,2),size(Y,1),1)-Y(:,2)).^2;
    
    mat=2*RBFscale^2*frbf(s,1)+4*RBFscale^4*syy.*frbf(s,2);   
    
    
 
    
        switch RBFtype
            case 'g'
   P=12*RBFscale^4*frbf(0,2)-lambda'*mat;
            case 'ms'
    P=12*RBFscale^4*frbf(0,2)-2*lambda'*mat+...
        lambda'*frbf(RBFscale^2*DistanceMatrixSquare(Y,Y),0)*lambda;                
                
                
        end     
        
    
    
    case 'L'
        
  

    switch RBFtype
        case {'ms','mq' }
    P=laplaceXYkermat(X,X) -2*laplacekermat(X,Y)*lambda+...
    lambda'*frbf(RBFscale^2*DistanceMatrixSquare(Y,Y),0)*lambda;
    %P=laplaceXYkermat(X,X) -laplacekermat(X,Y)*lambda;    
    

    
        case 'g'
    s=RBFscale^2*DistanceMatrixSquare(X,Y);      
    mat=2*RBFscale^2*(dim*frbf(s,1)+2*s.*frbf(s,2)); 
    P=dim*(dim+2)*4*RBFscale^4-mat*lambda; 
    end
            
            
end

P=sqrt(P);
    