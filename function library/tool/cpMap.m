function [cp,cpnormal,label] = cpMap(ftpts,normals,P,indices) 
% version 4.0
    n = size(P,1);
    cp = zeros(size(P));
    cpnormal = zeros(size(P));
    label = zeros(size(P,1),1);
    for k = 1:n
        % k
        ftpts_sub = ftpts(indices(k,:),:);
        normals_sub = normals(indices(k,:),:);

        cp(k,:) = ftpts_sub(1,:);
        cpnormal(k,:)  = normals_sub(1,:);


        nref = normals_sub(1,:);
        pk = P(k,:);

        dist = norm(cp(k,:)  - pk);
        angle = acos((normals_sub*nref')/(norm(normals_sub)*norm(nref)));
        index2 = angle < pi/2;
        ptsLoc = rotate(bsxfun(@minus,ftpts_sub(index2,:),pk),nref,1);
        % Fit local reconstruction polynomial
        deg = 3;
        a = ls_polyfit(ptsLoc(:,1:end-1),ptsLoc(:,end),deg);
        % Newton's method for finding the minimizer
        xmin = newton(a,ptsLoc(1,1:end-1)');
        cp2 = rotate([xmin ls_polyval(a,xmin)],nref,2) + pk;
        dist2 = norm(cp2  - pk);
        if dist2 <= dist
            cp(k,:) = cp2;
            cpnormal(k,:) = geomNormals(a,xmin,nref);
        else
            label(k) = 1;
        end       
    end
    label = logical(label);
end