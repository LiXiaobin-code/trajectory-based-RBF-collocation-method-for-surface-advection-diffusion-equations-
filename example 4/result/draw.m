clc;
clear;
close;
load('Torus_11600.mat');
[Theta, Phi] = cartesian2torus(p(:,1),p(:,2),p(:,3), 1, 1/3);
alphac = 2;   
betac = 2;    
u0 = abs(sin(alphac.*Theta)) + abs(cos(betac*Phi)) - 0.7*abs(sin(alphac*Theta).*cos(betac*Phi));


figure(1);
load('TCNtorus,mat.mat')
trisurf(t,p(:,1),p(:,2),p(:,3),u);
colormap(jet)
h = colorbar;
clim([0,max(u(:,1))])
axis off
box off
shading interp
axis equal
zlim([min(p(:,3)) ,max(p(:,3))])
set(gca, 'Position', [0.05 0.05 0.85 0.9]);  % Adjust plot area to leave space for colorbar
set(h, 'Position', [0.9 0.25 0.03 0.5]);  % Adjust colorbar position
print('-dpng', '-r400', 'CN_Torus'); % Save the figure as a PNG file with 400 DPI


figure(2);
 load('SLtorus,mat.mat')
trisurf(t,p(:,1),p(:,2),p(:,3),u);
colormap(jet)
h = colorbar;
clim([0,max(u(:,1))])
axis off
box off
shading interp
axis equal
zlim([min(p(:,3)) ,max(p(:,3))])

set(gca, 'Position', [0.05 0.05 0.85 0.9]);  % Adjust plot area to leave space for colorbar
set(h, 'Position', [0.9 0.25 0.03 0.5]);  % Adjust colorbar position
print('-dpng', '-r400', 'CN_Torus2'); % Save the figure as a PNG file with 400 DPI


figure(3);

 load('TCNtorus,mat.mat')
trisurf(t,p(:,1),p(:,2),p(:,3),abs(u-u0));
colormap(jet)
h = colorbar;
h.Ruler.Exponent = -2;

clim([0,0.05])
axis off
box off
shading interp
axis equal
zlim([min(p(:,3)) ,max(p(:,3))])

set(gca, 'Position', [0.05 0.05 0.85 0.9]);  % Adjust plot area to leave space for colorbar
set(h, 'Position', [0.9 0.25 0.03 0.5]);  % Adjust colorbar position

print('-dpng', '-r400', 'SL_Torus'); % Save the figure as a PNG file with 400 DPI


figure(4);

load('SLtorus,mat.mat')
trisurf(t,p(:,1),p(:,2),p(:,3),abs(u-u0));
colormap(jet)
h = colorbar;
h.Ruler.Exponent = -2;
clim([0,0.05])
axis off
box off
shading interp
axis equal
zlim([min(p(:,3)) ,max(p(:,3))])
set(gca, 'Position', [0.05 0.05 0.85 0.9]);  % Adjust plot area to leave space for colorbar
set(h, 'Position', [0.9 0.25 0.03 0.5]);  % Adjust colorbar position
print('-dpng', '-r400', 'SL_Torus2'); % Save the figure as a PNG file with 400 DPI



% %%
% figure(5);
% % 环面参数
% R = 1;      % 主半径
% r = 1/3;      % 管道半径
% theta = linspace(0, 2*pi, 800);
% phi = linspace(0, 2*pi, 800);
% 
% [Theta, Phi] = meshgrid(theta, phi);
% 
% 
% % 环面参数化
% X = (R + r*cos(Theta)) .* cos(Phi);
% Y = (R + r*cos(Theta)) .* sin(Phi);
% Z = r*sin(Theta);
% 
% % 可视化环面
% % 生成连续非光滑图案
% alpha = 2;   % 波数
% beta = 2;    % 交错因子
% f = abs(sin(alpha*Theta)) + abs(cos(beta*Phi)) - 0.7*abs(sin(alpha*Theta).*cos(beta*Phi));
% colormap(jet); % 设置色彩映射为 jet
% mesh(Theta, Phi,f);
% % 定义x和y的值
% x_values = 0:pi/2:2*pi; % 0, pi/2, pi, 3*pi/2, 2*pi
% y_values = pi/4:pi/2:2*pi; % 0, pi/2, pi, 3*pi/2, 2*pi
% 
% xlim([-0., 2*pi + 0.]); % 设置x轴范围
% ylim([-0., 2*pi + 0.]); % 设置y轴范围
% 
% % 绘制竖直线
% for k = 1:length(x_values)
%     line([x_values(k) x_values(k)], ylim, 'Color', [0.3010, 0.7450, 0.9330], 'LineStyle', '--', 'LineWidth', 2);
% end
% 
% % 绘制水平线
% for k = 1:length(y_values)
%     line(xlim, [y_values(k) y_values(k)], 'Color', [0.3010, 0.7450, 0.9330], 'LineStyle', '--', 'LineWidth', 2);
% end
% 
% % 设置坐标轴的刻度标签
% xticks(0:pi/2:2*pi);
% xticklabels({'0', '\pi/2', '\pi', '3\pi/2', '2\pi'});
% yticks(pi/4:pi/2:2*pi);
% yticklabels({ '\pi/4', '3\pi/4', '5\pi/4', '7\pi/4'});
% 
% % 添加标签和标题
% xlabel('\theta');
% ylabel('\phi');
% % title('Grid Lines at Specified Intervals on XOY Plane');
% set(gca, 'Position', [0.1,0.1,0.9,0.9]); 
% hold off;
% print('-dpng', '-r400', 'initial_Torus2_theta_phi'); % Save the figure as a PNG file with 400 DPI
% 
