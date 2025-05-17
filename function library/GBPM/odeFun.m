function T = odeFun(ctrs,n,grid,d,y,dx)

y = reshape(y,size(ctrs));
[~,~,~,~,T] = resamplingGBPM(ctrs,n,grid,d,y,dx);

T = T(:);