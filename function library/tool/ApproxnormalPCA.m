function [ Normal ] = ApproxnormalPCA( X, num_pt )
%Z input point cloud
%num_pt: size of local stencil for approximation of normal

%Normal: n_Z*3
Normal=zeros(length(X),3);
index  = knnsearch(X ,X ,'k',num_pt);
for i=1:length(X) 
X_loc = X(index(i,:)',:);
n=surface_pca(X_loc);
if n(end)~= 0
Normal(i,:) = sign(n(end))*n;
else
Normal(i,:) = n;  
end
%S=sqrt(diag(S));
%Curvature(i) = S(end)./sum(S);

end

