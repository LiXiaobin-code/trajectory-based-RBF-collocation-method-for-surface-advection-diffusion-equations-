function [ftpts,nref,flag] = collectPoints(ctrs,n,p,nref,m,options)
% This function finds a number of points on the surface based on a
% neighborhood of the reference grid point. 
%
% Input: 
%   activGPts        : The active Cartesian grid points.
%   ctrs             : The current points on the surface.
%   n                : The unit normal vector of the points on the surface.
%   p                : The reference grid point
%   m                : The number of points on the surface for the local reconstruction.
%   dx               : The grid size.
%
% Output:
%   ftpts    : The m closest points on the surface.
%   n0       : The unit normal vector of the closest points.
%   flag     : 1 if the function is able to find m points. 

dim = size(ctrs,2);

% The minimum distance of the points on the surface.
minDist = options.minDist;

% Sort the corresponding footpoints according to the distance from the
% reference grid point p.
[~,N2] = sort(sqrt(sum(bsxfun(@minus,ctrs,p).^2,2)));
ctrs = ctrs(N2,:); n = n(N2,:);

% Initialization
ftpts = zeros(m,dim);
k = 1;
ftpts(1,:) = ctrs(1,:); 
flag = 1;

% Loop for all the points in the neighborhood (x,y).
for l = 1:length(ctrs(:,1))
    % Find the minimum distance between the candidate point and the points
    % collected so far.
    if k>1
        dis = min(sqrt(sum(bsxfun(@minus,ftpts(1:k-1,:),ctrs(l,:)).^2,2)));
    else
        dis = 1e10;
    end
    
    % Find the angle between the candidate point and the reference normal
    % vector.
    angle = acos((n(l,:)*nref')/(norm(n(l,:))*norm(nref)));
    
    if angle > options.angleDeactiv
        flag = 0;
        break;
    end
    
    % Check for minimum distance and consistent Lagrangian information
    if dis>minDist && angle<options.maxAngle
        ftpts(k,:) = ctrs(l,:);
        k = k+1;
    end
    
    % Check whether m points are collected
    if k>m
        flag = 1;
        break;
    else
        flag = 0;
    end
end