function [ctrs,n,curv,d,tang1,tang2] = resamplingGBPM(ctrs,n,Xnew,dx,options)
% It returns the updated closest points of the surface for the given grid
% points.
%
% Input: 
%   ctrs    : The current points on the surface.
%   n       : The unit normal vector of the points on the surface.
%   Xnew    : The points to map on the surface.
%   dx      : The grid size.
%   options : struct that consists of the following:
%                   - modified : 1 if GBPM modified is used, 0 otherwise
%                   - deg      : polynomial degree for local reconstruction
%                                   (default value 2)
%
% Output:
%   ctrs   : The resampled points on the surface.
%   n      : The unit normal vector of the resampled surface points.
%   curv   : The curvature at the resampled surface points. 
%   d      : The active grid points after the resampling step.

% 它返回给定网格点的表面上最近点的更新结果。
%
% 输入:
%   ctrs    : 表面上当前的点。
%   n       : 表面点的单位法向量。
%   Xnew    : 要映射到表面的点。
%   dx      : 网格大小。
%   options : 结构体，包括以下内容：
%                   - modified : 如果使用修改版GBPM则为1，否则为0
%                   - deg      : 用于局部重建的多项式度数
%                                   (默认值为2)
%
% 输出:
%   ctrs   : 表面上重采样的点。
%   n      : 重采样表面点的单位法向量。
%   curv   : 重采样表面点的曲率。
%   d      : 重采样步骤后的活动网格点。


dim = size(n,2);
d2 = (1:length(Xnew(:,1)))';
if nargin < 5
    options.deg = 2;
    options.angleDeactiv = pi;
    options.maxAngle = pi/2;
    options.minDist = sqrt(dim)*dx/5;
end
deg = options.deg; % The degree of the local reconstruction polynomial.

if dim == 2
    P = deg + 1; % Number of polynomial terms
    m = round(2.5*P); % Number of points for surface reconstruction.
elseif dim == 3
    P = round((deg+1)*(deg+2)/2); % Number of polynomial terms
    m = round(2.5*P); % Number of points for surface reconstruction.
else
    error('Not implemented for higher dimensions')
end

% Initialize variables
newCtrs = zeros(length(d2),dim); newN = newCtrs;  tang1 = newCtrs; tang2 = newCtrs;
newCurv = zeros(length(d2),dim-1);

% Flag that points to active grid points after the resampling process.
active = zeros(size(d2));

% Loop for all the grid points for the resampling process
parfor j = 1:length(d2)
    p = Xnew(d2(j),:); % Reference grid point
    nref = n(d2(j),:); % Reference normal
    
    % Collect the no number of points for the local reconstruction. 
    [ftpts, nref, flag] = collectPoints(ctrs, n, p, nref, m, options);
    
    if flag % We successfully found points for local reconstruction.
        
        % Convert to local coordinates
        ptsLoc = rotate(bsxfun(@minus,ftpts,p),nref,1);
        
        % Fit local reconstruction polynomial
        a = ls_polyfit(ptsLoc(:,1:end-1),ptsLoc(:,end),deg);
        
        % Newton's method for finding the minimizer
        xmin = newton(a,ptsLoc(1,1:end-1)');
        
        % Calculate geometric quantities
        [nMin,curvMin,tangMin] = geomQuantities(a,xmin,nref);
        
        newCtrs(j,:) = rotate([xmin ls_polyval(a,xmin)],nref,2)+p;
        newN(j,:) = nMin;
        newCurv(j,:) = curvMin;
        tang1(j,:) = tangMin(1,:);
        if dim == 3
            tang2(j,:) = tangMin(2,:);
        end
        
        % Inclusion of points check
        if sum(sign((max(ptsLoc(:,1:end-1)) - xmin).*(xmin - min(ptsLoc(:,1:end-1)))))>length(xmin)-1
            active(j) = 1;
        end
    end
end
active = logical(active);
if sum(active) < length(active)
    disp('Points deactivated')
end
ctrs = newCtrs(active,:); n = newN(active,:); curv = newCurv(active,:);
d = d2(active);