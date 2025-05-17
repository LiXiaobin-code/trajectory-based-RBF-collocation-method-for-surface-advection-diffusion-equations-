clc;
clear;
close all;
% SL

global RBFtype RBFscale RBFpar;
RBFtype = 'ms';
RBFscale = 2;
dim = 2;
time_steps = [0.0005, 0.001,0.005,0.01 ]; % Example time step values

epsilon = 1e-6; % Fixed epsilon value

for dt_idx = 1:length(time_steps)
    dt = time_steps(dt_idx);
    
    if strcmp(RBFtype, 'ms')
        m = 4;
        RBFpar = m - dim / 2;

    end
    
    n = 500;
    ctrs = pickpointscircle(n);
    dsites = pickpointscircle(n);
    syms x y t;
    v = @(x)  [-x(:, 2), x(:, 1)];

    theta = linspace(0, 2 * pi, length(dsites) + 1)';
    theta(end) = [];
    u0 = initial(theta);

    Tf = 2*pi;
    nt = ceil(Tf / dt) + 1;
    nctrs = length(ctrs);
    alpha = zeros(nctrs, nt);

    [IM, surflap] = KPMDiffSurfMatrix(dsites, ctrs, 'Circle');
    CN.K = IM - dt / 2 * epsilon * surflap;

    v0 = v(dsites);
    Xback1 = dsites - dt / 2 * v0;
    ptsBack1 = cpCircle(Xback1);
    v1 = v(ptsBack1);
    Xback = dsites - dt * v1;
    ptsBack = cpCircle(Xback);

    [IM_Back, surflap_Back] = KPMDiffSurfMatrix(ptsBack, ctrs, 'Circle');
    CN.Kl1 = IM_Back + dt / 2 * epsilon * surflap_Back;

    [CN.Q, CN.R] = qr(CN.K);
    optsu.UT = true;

    alpha(:, 1) = IM\u0;
    u(:, 1) = u0;
    for i = 2:nt
        U = CN.Kl1 * alpha(:, i - 1);
        alpha(:, i) = linsolve(CN.R, CN.Q' * U, optsu);
        u(:, i) = IM * alpha(:, i);
        figure(1);
        clf;
        plot(theta, u(:, i), 'r.-');
    end

    SLCN = u(:, end);
    filename = sprintf('SLCN_dt_%g.mat', dt);
    save(filename, 'SLCN', 'theta', 'u', 'nt');
end
