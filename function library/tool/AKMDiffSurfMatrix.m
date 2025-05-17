function [Kmat,surflap]= AKMDiffSurfMatrix(X, Z,normal)
%version 3.0
    global RBFscale;
    [nx ,dim] = size(X);
    nz = size(Z,1);
    P = cell(3,3);
    for i=1:3
        for j=i:3
            P{i,j}=sparse(nx,nx); 
        end
    end
    for i=1:nx
        P{1,1}(i,i) = 1 - normal(i,1).*normal(i,1);
        P{1,2}(i,i) =   - normal(i,1).*normal(i,2);
        P{1,3}(i,i) =   - normal(i,1).*normal(i,3);
        P{2,2}(i,i) = 1 - normal(i,2).*normal(i,2);
        P{2,3}(i,i) =   - normal(i,2).*normal(i,3);
        P{3,3}(i,i) = 1 - normal(i,3).*normal(i,3);
    end

    Kmat = frbf(DistanceMatrixSquare(X,Z)*diag(RBFscale.^2),0);
    Euclgrad = zeros(nx,nz,dim);
    fmat1 = frbf(DistanceMatrixSquare(X,Z)*diag(RBFscale.^2),1)*diag(RBFscale.^2);
    surfgrad = zeros(nx,nz,dim);
    for i = 1:dim
        Euclgrad(:,:,i) = 2*fmat1.*(repmat(X(:,i),1,nz)-repmat(Z(:,i)',nx,1));
    end

    surfgrad(:,:,1) = (P{1,1})*Euclgrad(:,:,1) + (P{1,2})*Euclgrad(:,:,2) + (P{1,3})*Euclgrad(:,:,3);
    surfgrad(:,:,2) = (P{1,2})*Euclgrad(:,:,1) + (P{2,2})*Euclgrad(:,:,2) + (P{2,3})*Euclgrad(:,:,3);
    surfgrad(:,:,3) = (P{1,3})*Euclgrad(:,:,1) + (P{2,3})*Euclgrad(:,:,2) + (P{3,3})*Euclgrad(:,:,3);

    clear P fmat1 Euclgrad ;

    [Q,R] = qr(Kmat,"econ");
    opts.UT = true;
    temp = Q'*surfgrad(:,:,1);
    surflap = surfgrad(:,:,1)*linsolve(R, temp, opts);

    temp = Q'*surfgrad(:,:,2);
    surflap =  surflap +  surfgrad(:,:,2)*linsolve(R, temp, opts);

    temp = Q'*surfgrad(:,:,3);
    surflap = surflap + surfgrad(:,:,3)*linsolve(R, temp, opts);


    % surflap = surfgrad(:,:,1)/Kmat*surfgrad(:,:,1) + surfgrad(:,:,2)/Kmat*surfgrad(:,:,2) + surfgrad(:,:,3)/Kmat*surfgrad(:,:,3);

end


