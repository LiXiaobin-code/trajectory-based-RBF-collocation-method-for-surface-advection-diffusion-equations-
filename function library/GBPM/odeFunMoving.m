function T = odeFunMoving(ctrs,n,grid,d,y,dx,t)

y = reshape(y,size(ctrs));
[~,~,~,~,T] = resamplingGBPM(ctrs,n,grid,d,y,dx);

T = - T(:) - 0.05*cos(t);