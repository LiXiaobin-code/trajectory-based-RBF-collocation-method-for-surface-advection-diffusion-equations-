function normal_loc = mynormals(p,t,npts)

    % plot3(dsites(:,1),dsites(:,2),dsites(:,3),'r.')
    normals = vnormals(p, t);
    % Int_S = dsites - 1e-2*normals;
    % Ext_S = dsites + 1e-2*normals;
    % plot3(Int_S(:,1),Int_S(:,2),Int_S(:,3),'b.')
    % hold on;
    % patch('Vertices', dsites, 'Faces', faces, 'FaceVertexCData', 0*dsites(:,3), 'FaceColor', 'interp', 'EdgeColor', 'none');
    % plot3(Ext_S(:,1),Ext_S(:,2),Ext_S(:,3),'r.')

    Mdl = KDTreeSearcher(p);
    [indices, ~] = knnsearch(Mdl, p, 'K', npts); % 包括点本身
    ndsites = size(p,1);
    normal_loc = zeros(size(p));
    
    
    global RBFtype RBFscale RBFpar;
    
    RBFtype = 'ms';
    RBFscale = 1.;
    dim = 3;
    if strcmp(RBFtype,'ms')
        m = 4;
        RBFpar = m - dim/2;
    end
    
    
    for i = 1:ndsites
        i
        locdsites = p(indices(i,:),:);
        locnormals= normals(indices(i,:),:);
        locInt_S = locdsites - 1e-2*locnormals;
        locExt_S = locdsites + 1e-2*locnormals;
    
        S_interp = [locInt_S; locdsites ;locExt_S];
        S = S_interp;
        ns = size(S,1);
        b = zeros(size(S_interp,1),1); b(1:ns/3) = -1; b(2/3*ns+1:ns) = 1;
        IM = frbf(DistanceMatrixSquare(S_interp,S)*diag(RBFscale.^2),0);
        lam_surf = IM\b;
        X = p(i,:);
        nx = size(X,1);
        fmat1 = frbf(DistanceMatrixSquare(X,S)*diag(RBFscale.^2),1)*diag(RBFscale.^2);
        normal = zeros(nx,3);
        D1mat = zeros(nx,ns,dim);
        for j = 1:dim
            D1mat(:,:,j) = 2*fmat1.*(repmat(X(:,j),1,ns)-repmat(S(:,j)',nx,1));
            normal(:,j) = D1mat(:,:,j)*lam_surf;
        end
        normal_loc(i,:) = normalize(normal, 2, 'norm');
    end

end