clc;
clear;
close;

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

[IM, surflap, surfgrad]  =  KPMDiffSurfMatrix(dsites,ctrs,'Torus',2);
V = v(dsites);
G = surfgrad(:,:,1).*V(:,1) + surfgrad(:,:,2).*V(:,2) + surfgrad(:,:,3).*V(:,3);
%%
% CN
CN.K = IM - epsilon*dt/2*surflap;

[CN.Q,CN.R] = qr(CN.K);  
optsu.UT = true;
alpha(:,1) = IM\u0;

nt = ceil(Tf/dt)+1;

CN.K = IM + dt/2*G - dt/2*epsilon*surflap;
CN.Kl1 = IM - dt/2*G + dt/2*epsilon*surflap;
[CN.Q,CN.R] = qr(CN.K);  
optsu.UT = true;

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
save('TCNtorus,mat','u','alpha');
