clc;
clear;
close all;
% SL
global RBFtype RBFscale RBFpar;
RBFtype = 'ms';
RBFscale =  2;
dim = 2;
if strcmp(RBFtype, 'ms')
    m = 4;
    RBFpar = m - dim / 2;
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
dt = 1 /100;
nt = ceil(Tf / dt) + 1;
nctrs = length(ctrs);
alpha = zeros(nctrs, nt);

v0 = v(dsites);
Xback1 = dsites - dt / 2 * v0;
ptsBack1 = cpCircle(Xback1);
v1 = v(ptsBack1);
Xback = dsites - dt * v1;
ptsBack = cpCircle(Xback);

[IM, surflap] =  KPMDiffSurfMatrix(dsites, ctrs, 'Circle');
[IM_Back, surflap_Back] =  KPMDiffSurfMatrix(ptsBack, ctrs, 'Circle');



epsilons = [ 1e-3, 1e-4,1e-5, 1e-6]; % Example epsilon values

for eps_idx = 1:length(epsilons)
    epsilon = epsilons(eps_idx);

    CN.K = IM - dt / 2 * epsilon * surflap;
    [CN.Q, CN.R] = qr(CN.K);
    optsu.UT = true;

    CN.Kl1 = IM_Back + dt / 2 * epsilon * surflap_Back;
    alpha(:, 1) = IM \ u0;
    u(:, 1) = u0;
    for i = 2:nt
        U = CN.Kl1 * alpha(:, i - 1);
        alpha(:, i) = linsolve(CN.R, CN.Q' * U, optsu);
        u(:, i) = IM * alpha(:, i);
        % u(:, i)= imgaussfilt(u(:, i), 2);
        figure(1);
        clf;  
        plot(theta, u(:, i), 'r.-');
        title(sprintf('\\epsilon = %5.0e, T = %.2f',epsilon, i*dt));
        drawnow;
    end
    SLCN = u(:, end);
    filename = sprintf('SLCN_epsilon_%.0e.mat', epsilon);
    save(filename, 'SLCN');
end
