clc;
clear;
close;

global RBFtype RBFscale RBFpar;
RBFtype = 'ms';
RBFscale = 5;
dim = 2;
if strcmp(RBFtype,'ms')
    m = 6;
    RBFpar = m - dim/2;
end

syms x y t;
v = @(x)[-x(:,2),x(:,1)];


epsilon = 1e-1;
f = matlabFunction(-cos(x)*cos(y)*sin(t)-y*(-sin(x)*cos(y)*cos(t)/(x^2+y^2)*y^2+cos(x)*sin(y)*cos(t)/(x^2+y^2)*x*y)+x*(sin(x)*cos(y)*cos(t)/(x^2+y^2)*x*y-cos(x)*sin(y)*cos(t)/(x^2+y^2)*x^2)+epsilon*(((x^2+y^2)*cos(y)-y*sin(y))*cos(x)+2*x*(y*sin(y)-1/2*cos(y))*sin(x))*cos(t)/(x^2+y^2));

Tf = 1;
NIter = [50,65,80,100,120,150,200];
E = zeros(length(NIter),1);
for k = 1:length(NIter) 
    k
    dt = 1e-4;
    nt = floor(Tf/dt);

    n = NIter(k);
    ctrs = pickpointscircle(n);
    dsites = ctrs;
    uex = @(t) cos(dsites(:,1)).*cos(dsites(:,2)).*cos(t);
    u0 = uex(0);

    nctrs = length(ctrs);
    alpha = zeros(nctrs,nt);
    err = zeros(nt,1);

    [IM, surflap]  =  KPMDiffSurfMatrix(dsites,ctrs,'Circle');

    

    alpha(:,1) = IM\u0;


    SBDF3.K = 11/6*IM - dt*epsilon*surflap; 
    [SBDF3.Q,SBDF3.R] = qr(SBDF3.K,0);
    optsu.UT = true;

    %%
    % SBDF3

    SBDF3.K = 11/6*IM - dt*epsilon*surflap; 
    [SBDF3.Q,SBDF3.R] = qr(SBDF3.K,0);
    optsu.UT = true;
    %%
    v0 = v(dsites);
    Xback1 = dsites -  dt/2*(v0);
    ptsBack1 = cpCircle(Xback1);

    v1 = v(ptsBack1);
    Xback2 = dsites - dt*(-(v0) + 2*(v1));
    ptsBack2 = cpCircle(Xback2);

    v2 = v(ptsBack2);            
    Xback = dsites - dt/6*(v0 + 4*v1 + v2);
    ptsBack = cpCircle(Xback);
    pts_Back = ptsBack;
    SBDF3.Kl1 = frbf(DistanceMatrixSquare(ptsBack,ctrs)*diag(RBFscale.^2),0);
    [IM_Back, surflap_Back]  =  KPMDiffSurfMatrix(ptsBack,ctrs,'Circle');
    %%
    v0 = v(dsites);
    Xback1 = dsites -  2*dt/2*(v0);
    ptsBack1 = cpCircle(Xback1);

    v1 = v(ptsBack1);
    Xback2 = dsites - 2*dt*(-(v0) + 2*(v1));
    ptsBack2 = cpCircle(Xback2);

    v2 = v(ptsBack2);            
    Xback = dsites - 2*dt/6*(v0 + 4*v1 + v2);
    ptsBack = cpCircle(Xback);
    SBDF3.Kl2 = frbf(DistanceMatrixSquare(ptsBack,ctrs)*diag(RBFscale.^2),0);
    %%
    v0 = v(dsites);
    Xback1 = dsites -  3*dt/2*(v0);
    ptsBack1 = cpCircle(Xback1);

    v1 = v(ptsBack1);
    Xback2 = dsites - 3*dt*(-(v0) + 2*(v1));
    ptsBack2 = cpCircle(Xback2);

    v2 = v(ptsBack2);            
    Xback = dsites - 3*dt/6*(v0 + 4*v1 + v2);
    ptsBack = cpCircle(Xback);
    SBDF3.Kl3 = frbf(DistanceMatrixSquare(ptsBack,ctrs)*diag(RBFscale.^2),0);



    for i = 2:nt
        if i == 2 
            U =  (IM_Back + dt/2*epsilon*surflap_Back)*alpha(:,i-1) + 0.5*dt*f(dt*(i-1),dsites(:,1),dsites(:,2)) + 0.5*dt*f(dt*(i-2),pts_Back(:,1),pts_Back(:,2));
            alpha(:,i) = (IM - dt/2*epsilon*surflap)\U;
        elseif i==3
            U = 2*SBDF3.Kl1*alpha(:,i-1) - 0.5*SBDF3.Kl2*alpha(:,i-2) + dt*f(dt*(i-1),dsites(:,1),dsites(:,2));
            alpha(:,i) = (1.5*IM - dt*epsilon*surflap)\U;
        else
            U = 18/6*SBDF3.Kl1*alpha(:,i-1) - 9/6*SBDF3.Kl2*alpha(:,i-2)+ 2/6*SBDF3.Kl3*alpha(:,i-3) + dt*f(dt*(i-1),dsites(:,1),dsites(:,2));
            alpha(:,i) = linsolve(SBDF3.R,SBDF3.Q'*U,optsu);
        end
        u = IM*alpha(:,i);
        err(i) = max(abs(u - uex((i-1)*dt)))/max(abs(uex((i-1)*dt)));
    end
    E(k) = err(end);
end



loglog(2*pi./NIter,E,'*-');
name = sprintf('3SBDF%d.mat',m);

save(name,'E');

