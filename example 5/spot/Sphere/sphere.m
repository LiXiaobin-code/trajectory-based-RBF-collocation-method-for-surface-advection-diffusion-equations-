clc;
clear;
close;

options.method = 'rk2';

global RBFtype RBFscale RBFpar;
RBFtype = 'ms';
RBFscale = 2;
dim = 3;
if strcmp(RBFtype,'ms')
    m = 4;
    RBFpar = m - dim/2;
end



surftype = 'Sphere';
n = 4482;

filename = sprintf('%s_%d.mat', surftype, n);
matz = load(filename);
matx = matz;

ctrs = matz.p;
nz = size(ctrs,1);
dsites = matx.p;
nx = size(dsites,1);

FaceX = matx.t;

Tf = 600; 
dt = 0.5;
step = ceil(Tf/dt)+1;


epsil_w = 6.5e-3;
epsil_u = 0.516*epsil_w;
gam1 = 0.02; gam2 = 0.2;
a = 0.899; b = -0.91; c = -a;



[IM, surflap] = KPMDiffSurfMatrix(dsites,ctrs, surftype);

syms x y z;


switch surftype
    case 'Sphere' 
        R = 1;
        shape_f = @(x, y, z) x.^2 + y.^2 + z.^2 - R.^2;
    case 'Torus'
        R = 1;
        r = 1/3;
        shape_f = @(x, y, z) (x.^2 + y.^2 + z.^2 + R.^2 - r.^2).^2 - 4*R.^2.*(x.^2 + y.^2);
    case 'Bretzel2'
        shape_f = @(x, y, z) (-x.^4 + x.^2 - y.^2).^2 + z.^2/2 - (x.^2 + y.^2 + z.^2 + 1)/40;
    case 'Cyclide' 
        shape_f = @(x, y, z) (x.^2 + y.^2 + z.^2 - 1 + 19^2/100).^2 - 4*(2*x + sqrt(4 - 19^2/100)).^2 - (4*(19*y).^2)/100;
    case 'CPD surface'
        shape_f = @(x, y, z) sqrt((x - 1).^2 + y.^2 + z.^2).*sqrt((x + 1).^2 + y.^2 + z.^2).*sqrt((y - 1).^2 + x.^2 + z.^2).*sqrt((y + 1).^2 + x.^2 + z.^2) - 1.1;
    case 'orthocircle'
        shape_f = @(x, y, z) ((x.^2 + y.^2 - 1).^2 + z.^2).*((y.^2 + z.^2 - 1).^2 + x.^2).*((x.^2 + z.^2 - 1).^2 + y.^2) - 0.075.^2.*(1 + 3.*(x.^2 + y.^2 + z.^2));
end


N1 = diff(shape_f, x);
N2 = diff(shape_f, y);
N3 = diff(shape_f, z);
N{1} = matlabFunction(N1./sqrt(N1^2 + N2^2 + N3^2));
N{2} = matlabFunction(N2./sqrt(N1^2 + N2^2 + N3^2));
N{3} = matlabFunction(N3./sqrt(N1^2 + N2^2 + N3^2));

 
velocitytype = 1;

switch velocitytype
    case 1
        v = @(x)  0.01*[-N{2}(x(:,1),x(:,2),x(:,3)),N{1}(x(:,1),x(:,2),x(:,3)),0*x(:,3)];
    case 2
        v = @(x)  0.01*[0*N{1}(x(:,1),x(:,2),x(:,3)),-N{3}(x(:,1),x(:,2),x(:,3)),N{2}(x(:,1),x(:,2),x(:,3))];
    case 3
        vn = @(x)[N{1}(x(:,1),x(:,2),x(:,3)),N{2}(x(:,1),x(:,2),x(:,3)),N{3}(x(:,1),x(:,2),x(:,3))];
        vt1 = @(x) [x(:,1).*x(:,2),(x(:,1) + x(:,2)).*x(:,3),x(:,1).*x(:,2).*x(:,3)];
        v = @(x) 0.01*cross(vn(x),vt1(x),2);
end

switch options.method
    case 'rk1'
    %% RK1
        % Backward steps
        v0 = v(dsites); 
        Xback = dsites - dt*v0;
        ptsBack = cpSphere(Xback);
    case 'rk2'
    %% RK2
        % Backward steps
        v0 = v(dsites);  
        Xback1 = dsites -  dt/2*v0;
        ptsBack1 = cpSphere(Xback1);
        v1 = v(ptsBack1);
        Xback = dsites - dt*v1;
        ptsBack = cpSphere(Xback);
        SBDF1.Kl1 = frbf(DistanceMatrixSquare(ptsBack,ctrs)*diag(RBFscale.^2),0);

        v0 = v(dsites);
        Xback1 = dsites -  2*dt/2*v0;
        ptsBack1 = cpSphere(Xback1);     
        v1 = v(ptsBack1);
        Xback = dsites - 2*dt*v1;
        ptsBack = cpSphere(Xback);   
        SBDF2.Kl2 = frbf(DistanceMatrixSquare(ptsBack,ctrs)*diag(RBFscale.^2),0);
       
end


Au = 1.5*IM  - epsil_u*dt*surflap;
Aw = 1.5*IM  - epsil_w*dt*surflap;


uf = @(u,w) a*u.*(1-gam1*w.^2) + w.*(1-gam2*u);
wf = @(u,w) b*w.*(1 + a/b*gam1*u.*w) + u.*(c + gam2*w);

alpha = zeros(nz,step);
u =  zeros(nx,step);

beta = zeros(nz,step);
w =  zeros(nx,step);


stream1 = RandStream('mrg32k3a','seed',7122005);
uw = rand(stream1,nx,2)-0.5;
id = find(abs(dsites(:,3))>=0.1); 
uw(id,:) = 0; 
u0 = uw(:,1); 
w0 = uw(:,2);


alpha(:,1) = IM\u0;
u(:,1) = IM*alpha(:,1);
beta(:,1) = IM\w0;
w(:,1) = IM*beta(:,1);

[Qu,Ru] = qr(Au);
[Qw,Rw] = qr(Aw);
optsu.UT = true;

for i = 2:step
    i

    
    if i==2
        u1_back = SBDF1.Kl1*alpha(:,i-1);
        w1_back = SBDF1.Kl1*beta(:,i-1);

        rhsu = u1_back + dt*uf(u1_back,w1_back);
        rhsv = w1_back + dt*wf(u1_back,w1_back); 

        alpha(:,i) = (IM - epsil_u*dt*surflap)\rhsu;
        beta(:,i) = (IM - epsil_w*dt*surflap)\rhsv;  

        u(:,i) = IM*alpha(:,i);
        w(:,i) = IM*beta(:,i);
    else

        u1_back = SBDF1.Kl1*alpha(:,i-1);
        w1_back = SBDF1.Kl1*beta(:,i-1);

        u2_back = SBDF2.Kl2*alpha(:,i-2);
        w2_back = SBDF2.Kl2*beta(:,i-2);

        rhsu =  2*u1_back - 0.5*u2_back + dt*(2*uf(u1_back,w1_back) - uf(u2_back,w2_back));
        temp = Qu'*rhsu;
        alpha(:,i) = linsolve(Ru,temp,optsu);
           
        rhsv =  2*w1_back - 0.5*w2_back  + dt*(2*wf(u1_back,w1_back)-wf(u2_back,w2_back));
        temp = Qw'*rhsv;
        beta(:,i) = linsolve(Rw,temp,optsu);

        u(:,i) = IM*alpha(:,i); 
        w(:,i) = IM*beta(:,i);   
    end


    if(mod(i,10)==0)
        figure(1);
            trisurf(FaceX,dsites(:,1),dsites(:,2),dsites(:,3),u(:,i))
            axis equal;
            axis off;
            colormap(jet);
            view([-33, 55]);
            % 设置图形属性
            set(gca, 'Color', 'none');
            set(gcf, 'Color', 'white');  % 设置图形窗口的背景颜色为白色
            h.EdgeColor = 'none';
            h.LineStyle = 'none';
            shading interp;
            lighting phong;
            light('Position', [1 0 1], 'Style', 'infinite');
            h.FaceLighting = 'gouraud';
            h.AmbientStrength = 0.7;
            h.DiffuseStrength = 0.8;
            h.SpecularStrength = 0.9;
            h.SpecularExponent = 15;
            h.BackFaceLighting = 'unlit';
            set(gca, 'Unit', 'normalized', 'Position', [0, 0, 1, 1]);
            set(gcf, 'PaperPositionMode', 'auto');
            set(gcf, 'InvertHardcopy', 'off');
    end
end


load('Sphere_4482.mat')
filename = sprintf('u_%s_%d.mat', surftype, n);
save(filename,'u','p','t');






