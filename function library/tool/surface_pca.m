function [ normal ] = surface_pca(X)
%SURFACE_PCA Summary of this function goes here
%   Detailed explanation goes here
Dist  =X -mean(X ,1);
[~,~,V ]=svd(Dist);
normal=V (:,end);

end

