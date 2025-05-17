clc;
clear;
close;
%%
Tf = 6000; 
dt = 0.5;
step = ceil(Tf/dt)+1;


load('bunny15685.mat')
dsites = p;
faces = t;
npts = 100;
tic
normal_loc = mynormals(p,t,npts);
toc

ctrs = dsites;
normals = normal_loc;
Mdl = KDTreeSearcher(dsites);

nctrs = size(ctrs,1);
ndsites = size(dsites,1);

N = 20;
vscale = 0.01;
%% RK2
% Backward steps
v0 = vscale*cross(normals,repmat([0,0,1],ndsites,1),2);
Xback1 = dsites -  dt/2*v0;
[indices,~] = knnsearch(Mdl, Xback1, 'K', N); 
[~,cpnormal,label] = cpMap(dsites,normals,Xback1,indices);
v1 = vscale*cross(cpnormal, repmat([0,0,1],ndsites,1),2);  



Xback = dsites - dt*v1;
[indices,~] = knnsearch(Mdl, Xback, 'K', N); 
[ptsBack,ptsBacknormal] = cpMap(dsites,normals,Xback,indices);
v2 = vscale*cross(ptsBacknormal,repmat([0,0,1],ndsites,1),2);


Xback1 = ptsBack -  dt/2*v2;
[indices,~] = knnsearch(Mdl, Xback1, 'K', N); 
[~,cpnormal] = cpMap(dsites,normals,Xback1,indices);


v3 = vscale*cross(cpnormal, repmat([0,0,1],ndsites,1),2);
Xback2 = ptsBack - dt*v3;
[indices,~] = knnsearch(Mdl, Xback2, 'K', N); 
[ptsBack2,ptsBacknormal2] = cpMap(dsites,normals,Xback2,indices);

global RBFtype RBFscale RBFpar;

RBFtype = 'ms';
RBFscale = 2.;
dim = 3;
if strcmp(RBFtype,'ms')
    m = 4;
    RBFpar = m - dim/2;
end
[IM, surflap]= AKMDiffSurfMatrix(dsites,ctrs,normals);
SBDF2.Kl1 = frbf(DistanceMatrixSquare(ptsBack,ctrs)*diag(RBFscale.^2),0);  
SBDF2.Kl2 = frbf(DistanceMatrixSquare(ptsBack2,ctrs)*diag(RBFscale.^2),0);




% epsil_w = 4.5e-4;
epsil_w = 8.87e-4;
epsil_u = 0.516*epsil_w;
gam1 = 3; gam2 = 0;
a = 0.899; b = -0.91; c = -a;
uf = @(u,w) a*u.*(1-gam1*w.^2) + w.*(1-gam2*u);
wf = @(u,w) b*w.*(1 + a/b*gam1*u.*w) + u.*(c + gam2*w);


Au = 1.5*IM  - epsil_u*dt*surflap;
Aw = 1.5*IM  - epsil_w*dt*surflap;
alpha = zeros(nctrs,step);
% u =  zeros(ndsites,step);
beta = zeros(nctrs,step);
% w =  zeros(ndsites,step);
% init = @(x,y,z) atan(5*(-x-1.2));
% u(:,1) = init(dsites(:,1),dsites(:,2),dsites(:,3));
% w(:,1)= u(:,1);

stream1 = RandStream('mrg32k3a','seed',7122005);
uw = rand(stream1,ndsites,2)-0.5;
id = find(abs(dsites(:,3))>=0.1); 
uw(id,:) = 0; 
% u(:,1) = uw(:,1); 
% w(:,1) = uw(:,2);


alpha(:,1) = IM\uw(:,1);
beta(:,1) = IM\uw(:,2);
optsL.LT = true;  % L 是下三角矩阵
optsU.UT = true;  % U 是上三角矩阵
% [Lu, Uu, pu] = lu(Au,'vector');
% [Lw, Uw, pw] = lu(Aw,'vector');

[Qu,Ru] = qr(Au);
[Qw,Rw] = qr(Aw);
optsu.UT = true;


clear Au Aw;
for i = 2:step
    i
    if i==2
        u1_back = SBDF2.Kl1*alpha(:,i-1);
        w1_back = SBDF2.Kl1*beta(:,i-1);

        rhsu = u1_back + dt*uf(u1_back,w1_back);
        rhsw = w1_back + dt*wf(u1_back,w1_back); 

        alpha(:,i) = (IM - epsil_u*dt*surflap)\rhsu;
        beta(:,i) = (IM - epsil_w*dt*surflap)\rhsw;  

        % u(:,i) = IM*alpha(:,i);
        % w(:,i) = IM*beta(:,i);
    else
        u1_back = SBDF2.Kl1*alpha(:,i-1);
        w1_back = SBDF2.Kl1*beta(:,i-1);
        u2_back = SBDF2.Kl2*alpha(:,i-2);
        w2_back = SBDF2.Kl2*beta(:,i-2);

        rhsu =  2*u1_back - 0.5*u2_back + dt*(2*uf(u1_back,w1_back) - uf(u2_back,w2_back));
        rhsw =  2*w1_back - 0.5*w2_back  + dt*(2*wf(u1_back,w1_back)-wf(u2_back,w2_back));
        temp = Qu'*rhsu;
        alpha(:,i) = linsolve(Ru,temp,optsu);
        temp = Qw'*rhsw;
        beta(:,i) = linsolve(Rw,temp,optsu);

        % alpha(:,i) = linsolve(Uu, linsolve(Lu,  rhsu(pu), optsL), optsU);
        % beta(:,i) = linsolve(Uw, linsolve(Lw, rhsw(pw), optsL), optsU);
        % u(:,i) = IM*alpha(:,i); 
        % w(:,i) = IM*beta(:,i);   
    end
    if(mod(i,10)==0)
        figure(1);
        temp = IM*alpha(:,i);
        trisurf(faces,dsites(:,1),dsites(:,2),dsites(:,3),temp)
        colormap(jet)
        colorbar
        shading interp
        axis equal
        zlim([min(dsites(:,3)) - 0.1,max(dsites(:,3)) + 0.1])  
        title(sprintf('Time %d', i*dt));  
    end
end



surftype = 'bunny';
n = 15685;
filename = sprintf('u_%s_%d.mat', surftype, n);
u = IM*alpha;
save(filename,'u','p','t');