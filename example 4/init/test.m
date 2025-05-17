close all;
clc;
clear;
output_folder = 'C:\Users\Admin\Desktop\fig';

matx = load('Torus_11600.mat');
dsites = matx.p;
FaceX = matx.t;
R = 1;      
r = 1/3;  
alphac = 2;   
betac = 2; 


fig = figure(1);
[Theta, Phi] = cartesian2torus(dsites(:,1),dsites(:,2),dsites(:,3), R, r);   
f = abs(sin(alphac.*Theta)) + abs(cos(betac*Phi)) - 0.7*abs(sin(alphac*Theta).*cos(betac*Phi));
trisurf(FaceX,dsites(:,1),dsites(:,2),dsites(:,3),abs(f));
clim([0,1.3876])
colormap(jet)




shading interp

axis equal
axis off;
box off
zlim([min(dsites(:,3)) - 0.1,max(dsites(:,3)) + 0.1])

% print('-dpng', '-r400', 'initial_Torus'); % Save the figure as a PNG file with 400 DPI
frame_file = sprintf('%s/initial_Torus.png', output_folder);
exportgraphics(fig, frame_file, 'Resolution', 200);

%%
theta = linspace(0, 2*pi, 800);
phi = linspace(0, 2*pi, 800);
[Theta, Phi] = meshgrid(theta, phi);
X = (R + r*cos(Theta)) .* cos(Phi);
Y = (R + r*cos(Theta)) .* sin(Phi);
Z = r*sin(Theta);
f = abs(sin(alphac*Theta)) + abs(cos(betac*Phi)) - 0.7*abs(sin(alphac*Theta).*cos(betac*Phi));

x_values = 0:pi/2:2*pi; % 0, pi/2, pi, 3*pi/2, 2*pi
y_values = pi/4:pi/2:2*pi; % 0, pi/2, pi, 3*pi/2, 2*pi

fig= figure;
colormap(jet); % 设置色彩映射为 jet
mesh(Theta, Phi,f);
xlim([-0., 2*pi + 0.]); % 设置x轴范围
ylim([-0., 2*pi + 0.]); % 设置y轴范围

% 绘制竖直线
for k = 1:length(x_values)
    line([x_values(k) x_values(k)], ylim, 'Color', [0.3010, 0.7450, 0.9330], 'LineStyle', '--', 'LineWidth', 2);
end

% 绘制水平线
for k = 1:length(y_values)
    line(xlim, [y_values(k) y_values(k)], 'Color', [0.3010, 0.7450, 0.9330], 'LineStyle', '--', 'LineWidth', 2);
end

% 设置坐标轴的刻度标签
xticks(0:pi/2:2*pi);
xticklabels({'0', '\pi/2', '\pi', '3\pi/2', '2\pi'});
yticks(pi/4:pi/2:2*pi);
yticklabels({ '\pi/4', '3\pi/4', '5\pi/4', '7\pi/4'});

% 添加标签和标题
xlabel('\theta');
ylabel('\phi');
% title('Grid Lines at Specified Intervals on XOY Plane');
hold off;
% print('-dpng', '-r400', 'initial_Torus_theta_phi'); % Save the figure as a PNG file with 400 DPI


frame_file = sprintf('%s/initial_Torus2_theta_phi.png', output_folder);
exportgraphics(fig, frame_file, 'Resolution', 200);