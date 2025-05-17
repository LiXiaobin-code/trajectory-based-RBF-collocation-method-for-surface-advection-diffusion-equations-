clc;
clear;
close all;
% CN
global RBFtype RBFscale RBFpar;
RBFtype = 'ms';
RBFscale = 2;
dim = 2;
time_steps =[0.01,0.0005, 0.001, 0.005]; % Example time step values

epsilon = 1e-6; % Fixed epsilon value

for dt_idx = 1:length(time_steps)
    dt = time_steps(dt_idx);
    
    if strcmp(RBFtype, 'ms')
        m =4;
        RBFpar = m - dim / 2;
    else
        RBFpar = 3;
    end
    
    n = 500;
    ctrs = pickpointscircle(n);
    dsites = pickpointscircle(n);
    syms x y t;
    v = @(x) [-x(:, 2), x(:, 1)];

    theta = linspace(0, 2 * pi, length(dsites) + 1)';
    theta(end) = [];
    u0 = initial(theta);

    Tf = 2*pi;
    nt = ceil(Tf / dt) + 1;
    nctrs = length(ctrs);
    alpha = zeros(nctrs, nt);

    [IM, surflap, surfgrad] = KPMDiffSurfMatrix(dsites, ctrs, 'Circle', 'all');
    
    V = v(dsites);
    G = surfgrad(:,:,1).*V(:,1) + surfgrad(:,:,2).*V(:,2);

    % CN
    CN.K = IM + dt/2*G - dt/2*epsilon*surflap;
    CN.Kl1 = IM - dt/2*G + dt/2*epsilon*surflap;
    [CN.Q, CN.R] = qr(CN.K);  
    optsu.UT = true;
    
    alpha(:,1) = IM \ u0;
    u(:,1) = u0; 
    
    for i = 2:nt
        U = CN.Kl1 * alpha(:,i-1);
        alpha(:,i) = linsolve(CN.R, CN.Q' * U, optsu);
        u(:,i) = IM * alpha(:,i);
        figure(1);
        clf;
        plot(theta, u(:,i), '.-');
    end
    
    tCN = u(:, end);
    filename = sprintf('tCN_dt_%g.mat', dt);
    save(filename, 'tCN', 'theta', 'u', 'nt');
end
