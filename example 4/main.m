clc;
clear;
close;

options.method = 'rk2';

global RBFtype RBFscale RBFpar;
RBFtype = 'ms';
RBFscale = 3;
dim = 3;
if strcmp(RBFtype,'ms')
    m = 6;
    RBFpar = m - dim/2;
end

R = 1;
r = 1/3;
v = @(x) velocityTorus(x,R,r);

Tf = 2*pi; 
dt = Tf/1000;
epsilon = 1e-6;


matx = load('Torus_11600.mat');
dsites = matx.p;
FaceX = matx.t;
ctrs = dsites;


[Theta, Phi] = cartesian2torus(dsites(:,1),dsites(:,2),dsites(:,3), R, r);
alphac = 2;   
betac = 2;    
u0 = abs(sin(alphac.*Theta)) + abs(cos(betac*Phi)) - 0.7*abs(sin(alphac*Theta).*cos(betac*Phi));


switch options.method
    case 'rk1'
    %% RK1
        % Backward steps
        v0 = v(dsites);
        Xback = dsites - dt*v0;
        ptsBack = cpTorus(Xback);
    case 'rk2'
    %% RK2
        % Backward steps
        v0 = v(dsites);
        Xback1 = dsites -  dt/2*v0;
        ptsBack1 = cpTorus(Xback1);
        ptsBackMid = ptsBack1;

        v1 = v(ptsBack1);
        Xback = dsites - dt*v1;
        ptsBack = cpTorus(Xback);
    case 'rk3'
    %% RK2 for midpoint
        % Backward steps
        v0 = v(dsites);
        Xback1 = dsites - dt/4*v0;
        ptsBack1 = cpTorus(Xback1);

        v1 = v(ptsBack1);
        Xback = dsites -  dt/2*v1;
        ptsBackMid = cpTorus(Xback);
    %% RK3
        % Backward steps
        v0 = v(dsites);
        Xback1 = dsites -  dt/2*v0;
        ptsBack1 = cpTorus(Xback1);

        v1 = v(ptsBack1);
        Xback2 = dsites - dt*(-v0 + 2*v1);
        ptsBack2 = cpTorus(Xback2);

        v2 = v(ptsBack2);            
        Xback = dsites - dt/6*(v0 + 4*v1 + v2);
        ptsBack = cpTorus(Xback);

    case 'rk4'
    %% RK3 for midpoint
        % Backward steps
        v0 = v(dsites);
        Xback1 = dsites - dt/4*v0;
        ptsBack1 = cpTorus(Xback1);

        v1 = v(ptsBack1);
        Xback2 = dsites -  dt/2*(-v0 + 2*v1);
        ptsBack2 = cpTorus(Xback2);

        v2 = v(ptsBack2);            
        Xback = dsites - dt/12*(v0 + 4*v1 + v2);
        ptsBackMid = cpTorus(Xback);
    %% RK4
        % Backward steps
        v0 = v(dsites);
        Xback1 = dsites -  dt/2*v0;
        ptsBack1 = cpTorus(Xback1);

        v1 = v(ptsBack1);            
        Xback2 = dsites -  dt/2*v1;
        ptsBack2 = cpTorus(Xback2);

        v2 = v(ptsBack2);            
        Xback3 = dsites - dt*v2;
        ptsBack3 = cpTorus(Xback3);

        v3 = v(ptsBack3);            
        Xback = dsites - dt/6*(v0 + 2*v1 + 2*v2 + v3);
        ptsBack = cpTorus(Xback);
end

[IM, surflap]  =  KPMDiffSurfMatrix(dsites,ctrs,'Torus');
[IM_Back, surflap_Back]  =  KPMDiffSurfMatrix(ptsBack,ctrs,'Torus');
%%
% CN
CN.K = IM - epsilon*dt/2*surflap;
CN.Kl1 = IM_Back + epsilon*dt/2*surflap_Back;


[CN.Q,CN.R] = qr(CN.K);  
optsu.UT = true;
alpha(:,1) = IM\u0;

nt = ceil(Tf/dt)+1;
for i = 2:nt
    i
    U = CN.Kl1*alpha(:,i-1);
    alpha(:,i) = linsolve(CN.R,CN.Q'*U,optsu);
    u = IM*alpha(:,i);
    if(mod(i,10)==0)
        figure(1);
        trisurf(FaceX,dsites(:,1),dsites(:,2),dsites(:,3),u)
        colormap(jet)
        colorbar
        shading interp
        axis equal
        zlim([min(dsites(:,3)) - 0.1,max(dsites(:,3)) + 0.1])  
    end

end
save('SLtorus,mat','u','alpha');
