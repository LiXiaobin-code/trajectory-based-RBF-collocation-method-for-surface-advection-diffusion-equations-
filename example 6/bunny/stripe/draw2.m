clc;
clear;
close;
load('u_bunny_15685.mat')

for i = 1:size(u,2)
    i
    if(mod(i,100)==0)

        figure(1);
        h = trisurf(t, p(:,1), p(:,2), p(:,3), u(:,i));
        axis equal;
        axis off;
        colormap(jet);
        % view([-33, 55]);
        
        % 设置图形属性
        set(gca, 'Position', [0 0 1 1]);
        set(gca, 'Color', 'none');
         set(gcf, 'Color', 'white');  % 设置图形窗口的背景颜色为白色
        h.EdgeColor = 'none';
        h.LineStyle = 'none';
        shading interp;
        lighting phong;
        light('Position', [1 -1 -5], 'Style', 'local');
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





% 创建图形并绘制三维曲面
% figure;
% h = trisurf(t, p(:,1), p(:,2), p(:,3), u(:,end));
% axis equal;
% axis off;
% colormap(jet);
% % view([-33, 55]);
% 
% % 设置图形属性
% set(gca, 'Color', 'none');
%  set(gcf, 'Color', 'white');  % 设置图形窗口的背景颜色为白色
% h.EdgeColor = 'none';
% h.LineStyle = 'none';
% shading interp;
% lighting phong;
% light('Position', [1 0 1], 'Style', 'infinite');
% h.FaceLighting = 'gouraud';
% h.AmbientStrength = 0.7;
% h.DiffuseStrength = 0.8;
% h.SpecularStrength = 0.9;
% h.SpecularExponent = 15;
% h.BackFaceLighting = 'unlit';
% set(gca, 'Unit', 'normalized', 'Position', [0, 0, 1, 1]);
% set(gcf, 'PaperPositionMode', 'auto');
% set(gcf, 'InvertHardcopy', 'off');
% 
% % 保存图像
% filename = sprintf('spot_%s_%d.png', surftype, n);
% print('-dpng', '-r400', filename);


fig = figure('Units', 'pixels', 'Position', [100 100 900 600]);
% output_folder = 'strip_sphere';
% if ~exist(output_folder, 'dir')
%     mkdir(output_folder);
% end
surftype = 'bunny';
n = 15685;
output_folder = 'C:\Users\Admin\Desktop\fig';
clf(fig);
h = trisurf(t, p(:,1), p(:,2), p(:,3), u(:,12001));
axis equal;
axis off;
colormap(jet);
% view([-33, 55]);

% % 添加光源
% light('Position', [1 -1 -5], 'Style', 'local');
% % 设置光照模型
% lighting gouraud;



% 设置图形属性
[cmin,cmax] = clim;
% set(gca,'CLim',[cmin + 0.025,cmax - 0.025]);
set(gca, 'Position', [0 0 1 1]);
set(gca, 'Color', 'none');
set(gcf, 'Color', 'white');  % 设置图形窗口的背景颜色为白色
h.EdgeColor = 'none';
h.LineStyle = 'none';
shading interp;
lighting phong;
% light('Position', [1 0 1], 'Style', 'infinite');
light('Position', [1 -1 -5], 'Style', 'local');

h.AmbientStrength = 0.7;
h.DiffuseStrength = 0.8;
h.SpecularStrength = 0.9;
h.SpecularExponent = 15;
h.BackFaceLighting = 'unlit';
set(gca, 'Unit', 'normalized', 'Position', [0, 0, 1, 1]);
set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'InvertHardcopy', 'off');
drawnow;
% pause(0.001);
frame_file = sprintf('%s/stripe_%s_%d.png', output_folder, surftype, n);
exportgraphics(fig, frame_file, 'Resolution', 400);














