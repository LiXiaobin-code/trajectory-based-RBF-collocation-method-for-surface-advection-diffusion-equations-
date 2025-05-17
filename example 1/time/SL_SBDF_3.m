clc;
clear;
close;

global RBFtype RBFscale RBFpar;
RBFtype = 'ms';
RBFscale = 2;
dim = 2;
if strcmp(RBFtype,'ms')
    m = 5;
    RBFpar = m - dim/2;
end
n = 200;
ctrs = pickpointscircle(n);
dsites = ctrs;
syms x y t;
v = @(x)[-x(:,2),x(:,1)];


epsilon = 1e-5;
f = matlabFunction(-cos(x)*cos(y)*sin(t)-y*(-sin(x)*cos(y)*cos(t)/(x^2+y^2)*y^2+cos(x)*sin(y)*cos(t)/(x^2+y^2)*x*y)+x*(sin(x)*cos(y)*cos(t)/(x^2+y^2)*x*y-cos(x)*sin(y)*cos(t)/(x^2+y^2)*x^2)+epsilon*(((x^2+y^2)*cos(y)-y*sin(y))*cos(x)+2*x*(y*sin(y)-1/2*cos(y))*sin(x))*cos(t)/(x^2+y^2));


uex = @(t) cos(dsites(:,1)).*cos(dsites(:,2)).*cos(t);


u0 = uex(0);
Tf = 1;
dtIter = [0.1,0.05,0.025,0.01,0.001];
E = zeros(length(dtIter),1);
for k = 1:length(dtIter) 
    dt = dtIter(k);
    nt = floor(Tf/dt);
    nctrs = length(ctrs);
    alpha = zeros(nctrs,nt);
    err = zeros(nt,1);

    [IM, surflap]  =  KPMDiffSurfMatrix(dsites,ctrs,'Circle');

    

    alpha(:,1) = IM\u0;


    SBDF3.K = 11/6*IM - epsilon*dt*surflap; 
    [SBDF3.Q,SBDF3.R] = qr(SBDF3.K,0);
    optsu.UT = true;

    %%
    % SBDF3

    SBDF3.K = 11/6*IM - epsilon*dt*surflap; 
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
            % alpha(:,i) = IM\uex(dt);
            % U = SBDF3.Kl1*alpha(:,i-1) + dt*f(dt*(i-1),dsites(:,1),dsites(:,2));
            % alpha(:,i) = (IM - dt*surflap)\U;
            U =  (IM_Back + epsilon*dt/2*surflap_Back)*alpha(:,i-1) + 0.5*dt*f(dt*(i-1),dsites(:,1),dsites(:,2)) + 0.5*dt*f(dt*(i-2),pts_Back(:,1),pts_Back(:,2));
            alpha(:,i) = (IM - epsilon*dt/2*surflap)\U;
        elseif i==3
            U = 2*SBDF3.Kl1*alpha(:,i-1) - 0.5*SBDF3.Kl2*alpha(:,i-2) + dt*f(dt*(i-1),dsites(:,1),dsites(:,2));
            alpha(:,i) = (1.5*IM - epsilon*dt*surflap)\U;
        else
            U = 18/6*SBDF3.Kl1*alpha(:,i-1) - 9/6*SBDF3.Kl2*alpha(:,i-2)+ 2/6*SBDF3.Kl3*alpha(:,i-3) + dt*f(dt*(i-1),dsites(:,1),dsites(:,2));
            alpha(:,i) = linsolve(SBDF3.R,SBDF3.Q'*U,optsu);
        end
        u = IM*alpha(:,i);
        err(i) = max(abs(u - uex((i-1)*dt)))/max(abs(uex((i-1)*dt)));
    end
    E(k) = err(end);
end


loglog(dtIter,E,'*-');


save('SBDF3.mat','E');

