function [ctrs,n,curv,d,tang1,tang2] = updateCompTubeGBPM(ctrs,n,tang1,tang2,curv,grid,M,d,bw,dx,options)

dim = size(n,2);

% Create a matrix corresponding to the Cartesian grid. The active grid
% points in the computational tube are labelled with the entry 2.
newBand = zeros(M);
newBand(d) = 2;

% Finds and adds all the grid points neighboring to active grid points and
% label them with the entry 1.
if dim == 2
    for L = 1:length(d)
        [k,l] = ind2sub(M,d(L));
        for p = k-1:k+1
            for q = l-1:l+1
                if p > 0 && p <= M(1) && q > 0 && q <= M(2)
                    if newBand(p,q)~=2 % Exclude points already in the band
                        newBand(p,q) = 1;
                    end   
                end
            end
        end
    end
elseif dim == 3
    for L = 1:length(d)
        [k,l,m] = ind2sub(M,d(L));
        for p = k-1:k+1
            for q = l-1:l+1
                for r = m-1:m+1
                    if p > 0 && p <= M(1) && q > 0 && q <= M(2) && r > 0 && r <= M(3)
                        if newBand(p,q,r)~=2 % Exclude points already in the band
                            newBand(p,q,r) = 1;
                        end
                    end
                end
            end
        end
    end
else
    error('Not implemented for higher dimensions')
end

% The newly activated grid points.
dWide = find(newBand(:)==1);
clear newBand

% A resampling step for the newly activated grid points.
[ctrsWide,nWide,curvWide,dWide,tang1Wide,tang2Wide] = resamplingGBPM(ctrs,n,grid,d,grid(dWide,:),dx,options);

% Add the new footpoints to the current ones.
ctrs = [ctrs;ctrsWide];
n = [n;nWide];
curv = [curv;curvWide];
d = [d;dWide];
tang1 = [tang1;tang1Wide];
tang2 = [tang2;tang2Wide];

% Find and deactivate the grid points that are far from the surface.
L = find(sqrt(sum((ctrs-grid(d,:)).^2,2))<=bw*dx); % Remove remote points
    
ctrs = ctrs(L,:);
n = n(L,:);
curv = curv(L,:);
d = d(L);
tang1 = tang1(L,:);
tang2 = tang2(L,:);