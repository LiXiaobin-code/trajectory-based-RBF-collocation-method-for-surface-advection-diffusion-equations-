clc;
clear;
close;

dtIter = [0.1,0.05,0.025,0.01,0.001];
A(:,1) = dtIter';
figure;
load('SBDF1.mat')
loglog(dtIter,E,'*-');
A(:,2) = E;
hold on;
load('CN.mat')
loglog(dtIter,E,'*-');
A(:,3) = E;
load('SBDF3.mat')
loglog(dtIter,E,'*-');
A(:,4) = E;
grid on;

legend('SL Back Euler','SL CN','SL SBDF3','Location','southeast');
% 创建 text
text('FontSize',16,'String','slope 1.0',...
    'Position',[0.0187830515937346 0.0804736915124337 0]);

% 创建 text
text('FontSize',16,'String','slope 2.0',...
    'Position',[0.0194198176253206 0.000664955371384009 0]);

% 创建 text
text('FontSize',16,'String','slope 3.0',...
    'Position',[0.0195581432597718 3.03543842909318e-06 0]);

xlabel('\Delta t');
ylabel('\infty-norm error')
print('a01.png', '-dpng', '-r300');

B = convergencerate(A);
mypath = fullfile(pwd,'timeconvergence.txt');
writertxt(B,mypath)


