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

    [IM, surflap]  = KPMDiffSurfMatrix(dsites,ctrs,'Circle');


        
   
    %%
    % SBDF1
    SBDF1.K = IM - epsilon*dt*surflap; 
    [SBDF1.Q,SBDF1.R] = qr(SBDF1.K,0);  
    optsu.UT = true;
    v0 = v(dsites);
    Xback1 = dsites -  dt*v0;
    ptsBack = cpCircle(Xback1);       
    SBDF1.Kl1 = frbf(DistanceMatrixSquare(ptsBack,ctrs)*diag(RBFscale.^2),0);

    alpha(:,1) = IM\u0;
    
    for i = 2:nt

            % i
            U = SBDF1.Kl1*alpha(:,i-1) + dt*f(dt*(i-1),dsites(:,1),dsites(:,2));
            alpha(:,i) = linsolve(SBDF1.R,SBDF1.Q'*U,optsu);
            u = IM*alpha(:,i);
            err(i) = max(abs(u - uex((i-1)*dt)))/max(abs(uex((i-1)*dt)));
            
    end
    
   
    E(k) = err(end);
end


loglog(dtIter,E,'*-');
save('SBDF1.mat','E');
