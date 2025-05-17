function [ctrs,n,curv,d,tang1,tang2] = GBPM(ctrs,n,grid,M,d,Xnew,V,T,dt,bw,dx,options)
% This function implements one full step of the Grid Based Particle 
% Method. 

if nargin < 12
    options.modified = 0;
    options.deg = 2;
    options.angleDeactiv = pi;
    options.maxAngle = pi/2;
end

[ctrs,n,~,d,~,~] = resamplingGBPM(ctrs,n,grid,d,Xnew,dx,options);

ctrs = motion(ctrs,V,n,T,dt);

[ctrs,n,curv,d,tang1,tang2] = resamplingGBPM(ctrs,n,grid,d,grid(d,:),dx,options);

[ctrs,n,curv,d,tang1,tang2] = updateCompTubeGBPM(ctrs,n,tang1,tang2,curv,grid,M,d,bw,dx,options);
