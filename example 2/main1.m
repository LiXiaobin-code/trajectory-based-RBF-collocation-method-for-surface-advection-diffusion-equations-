% clear;
% close all;
% 
% %% 设置采样间隔 k
% k = 1; % 恢复默认（不采样），若需采样可改为 k=10 等
% 
% %% Generate angular coordinates
% theta = linspace(0, 2*pi, 501)'; 
% theta(end) = [];
% theta_sampled = theta(1:k:end); % 采样后的 theta
% 
% %% Calculate initial condition
% u0 = initial(theta);
% u0_sampled = u0(1:k:end); % 采样后的 u0
% 
% epsilons = [1e-3, 1e-4, 1e-5, 1e-6];
% 
% % 自定义颜色（RGB格式，归一化到[0,1]）
% colors = [255, 187,  0;
%           246, 83,  20; 
%           124, 187,  0;
%           0, 161, 241] / 255;
% 
% % 线型组合：实线、虚线、点线、点划线
% line_styles = {'-', '--', ':', '-.'};
% 
% % 标记样式（可选，若需稀疏标记）
%  markers = {'none', 'none', 'none', 'none'}; % 无标记
%  % markers = {'o', 's', '^', 'd'}; % 圆形、方形、上三角、菱形
% 
% figure;
% set(gcf, 'Color', 'w'); % 白色背景
% hold on;
% 
% % 绘制初始条件（黑色圆圈）
%  % 绘制初始条件（黑色圆圈 + 填充色）
% plot(theta_sampled, u0_sampled, 'o', ...
%     'Color', [0 0 0], ...         % 边框颜色（纯黑）
%     'MarkerFaceColor', [0 0 0], ... % 填充色（白色）
%     'LineWidth', 1, ...
%     'MarkerSize', 2, ...
%     'DisplayName', 'Initial');
% 
% % 绘制不同 epsilon 的曲线
% for eps_idx = 1:length(epsilons)
%     epsilon = epsilons(eps_idx);
%     filename = sprintf('SLCN_epsilon_%.0e', epsilon);
%     load(filename); % 加载 SLCN
% 
%     SLCN_sampled = SLCN(1:k:end); % 采样后的 SLCN
% 
%     plot(theta_sampled, SLCN_sampled, ...
%         'Color', colors(eps_idx, :), ...
%         'LineStyle', line_styles{eps_idx}, ...
%         'Marker', markers{eps_idx}, ...
%         'MarkerSize', 3, ...
%         'LineWidth', 2, ...
%         'DisplayName', sprintf('\\epsilon = %.0e', epsilon));
% end
% 
% % 图例设置
% lgd = legend('show');
% lgd.Box = 'off';
% lgd.FontSize = 12;
% lgd.Location = 'best'; % 自动选择最佳位置
% 
% % 坐标轴标签
% xlabel(' \theta', 'FontSize', 16);
% ylabel('u', 'FontSize', 16);
% 
% % 坐标范围与刻度
% xlim([0, 2.*pi]);
% ylim([-0.15, 1.15]);
% set(gca, 'XTick', 0:pi/2:2*pi, ...
%          'XTickLabel', {'0', '\pi/2', '\pi', '3\pi/2', '2\pi'}, ...
%          'FontSize', 12);
% % grid on; % 添加网格线
% 
% % 可选：保存图像（600 dpi）
% print(gcf, '-dpng', '-r300', 'SLCN_epsilon_lines.png');

clear;
close all;

%% 设置采样间隔 k
k = 1; % 恢复默认（不采样），若需采样可改为 k=10 等

%% Generate angular coordinates
theta = linspace(0, 2*pi, 501)'; 
theta(end) = [];
theta_sampled = theta(1:k:end); % 采样后的 theta

%% Calculate initial condition
u0 = initial(theta);
u0_sampled = u0(1:k:end); % 采样后的 u0

epsilons = [1e-3, 1e-4, 1e-5, 1e-6];

% 自定义颜色（RGB格式，归一化到[0,1]）
colors = [255, 187,  0;
          246, 83,  20; 
          124, 187,  0;
          0, 161, 241] / 255;

% 线型组合：实线、虚线、点线、点划线
line_styles = {'-', '--', ':', '-.'};

% 标记样式（可选，若需稀疏标记）
markers = {'none', 'none', 'none', 'none'}; % 无标记
% markers = {'o', 's', '^', 'd'}; % 圆形、方形、上三角、菱形

%% 设置图像长宽比为1:2
figure;
set(gcf, 'Color', 'w', 'Units', 'inches', 'Position', [1 1 12 3]); % 宽度6英寸，高度3英寸（1:2比例）
hold on;

% 绘制初始条件（黑色圆圈 + 填充色）
plot(theta_sampled, u0_sampled, 'o', ...
    'Color', [0 0 0], ...         % 边框颜色（纯黑）
    'MarkerFaceColor', [0 0 0], ... % 填充色（黑色）
    'LineWidth', 1, ...
    'MarkerSize', 2, ...
    'DisplayName', 'Initial');

% 绘制不同 epsilon 的曲线
for eps_idx = 1:length(epsilons)
    epsilon = epsilons(eps_idx);
    filename = sprintf('SLCN_epsilon_%.0e', epsilon);
    load(filename); % 加载 SLCN
    
    SLCN_sampled = SLCN(1:k:end); % 采样后的 SLCN
    
    plot(theta_sampled, SLCN_sampled, ...
        'Color', colors(eps_idx, :), ...
        'LineStyle', line_styles{eps_idx}, ...
        'Marker', markers{eps_idx}, ...
        'MarkerSize', 3, ...
        'LineWidth', 2, ...
        'DisplayName', sprintf('\\epsilon = %.0e', epsilon));
end

% 图例设置：右上角
lgd = legend('show');
lgd.Box = 'off';
lgd.FontSize = 12;
lgd.Location = 'northeast'; % 将图例固定在右上角

% 坐标轴标签
xlabel(' \theta', 'FontSize', 16);
ylabel('u', 'FontSize', 16);

% 坐标范围与刻度
xlim([0, 2.*pi]);
ylim([-0.15, 1.15]);
set(gca, 'XTick', 0:pi/2:2*pi, ...
         'XTickLabel', {'0', '\pi/2', '\pi', '3\pi/2', '2\pi'}, ...
         'FontSize', 12);
% grid on
print(gcf, '-dpng', '-r600', 'SLCN_epsilon_lines.png');