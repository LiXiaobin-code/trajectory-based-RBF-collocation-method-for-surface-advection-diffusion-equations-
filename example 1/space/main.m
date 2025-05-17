clc;
clear;
close;
NIter = [50,65,80,100,120,150,200];


NIter =2*pi./NIter;
A(:,1) = NIter';
figure;
load('3SBDF3.mat')

loglog(NIter,E,'*-');A(:,2) = E;
hold on;
load('3SBDF4.mat')

loglog(NIter,E,'*-');A(:,3) = E;
load('3SBDF5.mat')

loglog(NIter,E,'*-');A(:,4) = E;
load('3SBDF6.mat')

loglog(NIter,E,'*-');A(:,5) = E;
grid on;

legend('m = 3','       4','       5','       6','Location','southeast');

xticks([0.04,0.06,0.08,0.1,0.12]);

xticklabels({'0.04','0.06','0.08','0.1','0.12'});


xlabel('\Delta x');
ylabel('\infty-norm error')

% 创建 text
text('FontSize',14,'String','slope 2.97',...
    'Position',[0.07 0.0138005493060751 0]);

% 创建 text
text('FontSize',14,'String','slope 4.98',...
    'Position',[0.07 3.99058227699502e-05 0]);

% 创建 text
text('FontSize',14,'String','slope 6.98',...
    'Position',[0.07 5.70790431199402e-07 0]);

% 创建 text
text('FontSize',14,'String','slope 8.95',...
    'Position',[0.07 3.06935154058372e-09 0]);

print('a02.png', '-dpng', '-r300');


B = convergencerate(A);
mypath = fullfile(pwd,'spatialconvergence.txt');
writertxt(B,mypath)