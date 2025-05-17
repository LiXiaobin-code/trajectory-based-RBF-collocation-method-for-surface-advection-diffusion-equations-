clc;
clear;
close;
 load('u_Sphere_4482.mat');

% for i = 1:size(u,2)
%     i
%     if(mod(i,10)==0)
%         % figure(1);
%         % trisurf(t,p(:,1),p(:,2),p(:,3),u(:,i))
%         % axis equal;
%         % axis off;
%         % colormap(jet);
%         % [cmin,cmax] = clim;
%         % climit = get(gca, 'CLim');
%         % % set(gca,'CLim',[cmin + 0.05,cmax - 0.05]);
%         % set(gca, 'Position', [0 0 1 1]);  % 将轴位置设置为填满整个图形窗口
%         % set(gca,'color','none')
%         % h.EdgeColor = 'none';
%         % h.LineStyle = 'none';
%         % % h.FaceColor=[1,0,0];
%         % shading interp ;
%         % lighting phong;
%         % light
%         % view([-33,55]);
%         % h.FaceLighting = 'gouraud';
%         % h.AmbientStrength = 0.7;
%         % h.DiffuseStrength = 0.8;
%         % h.SpecularStrength = 0.9;
%         % h.SpecularExponent = 15;
%         % h.BackFaceLighting = 'unlit';
% 
%         figure(1);
%         h = trisurf(t, p(:,1), p(:,2), p(:,3), u(:,i));
%         axis equal;
%         axis off;
%         colormap(jet);
%         view([-33, 55]);
% 
%         % 设置图形属性
%         set(gca, 'Position', [0 0 1 1]);
%         set(gca, 'Color', 'none');
%          set(gcf, 'Color', 'white');  % 设置图形窗口的背景颜色为白色
%         h.EdgeColor = 'none';
%         h.LineStyle = 'none';
%         shading interp;
%         lighting phong;
%         light('Position', [1 0 1], 'Style', 'infinite');
%         h.FaceLighting = 'gouraud';
%         h.AmbientStrength = 0.7;
%         h.DiffuseStrength = 0.8;
%         h.SpecularStrength = 0.9;
%         h.SpecularExponent = 15;
%         h.BackFaceLighting = 'unlit';
%         set(gca, 'Unit', 'normalized', 'Position', [0, 0, 1, 1]);
%         set(gcf, 'PaperPositionMode', 'auto');
%         set(gcf, 'InvertHardcopy', 'off');
% 
% 
% 
% 
% 
% 
% 
% 
%     end
% end

surftype = 'Sphere';
n = 4482;



% 创建图形并绘制三维曲面
fig = figure;
h = trisurf(t, p(:,1), p(:,2), p(:,3), u(:,end));
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

% % 保存图像
% filename = sprintf('spot_%s_%d.png', surftype, n);
% print('-dpng', '-r400', filename);

output_folder = 'C:\Users\Admin\Desktop\fig';
frame_file = sprintf('%s/spot_%s_%d.png', output_folder, surftype,n);
exportgraphics(fig, frame_file, 'Resolution', 200);


