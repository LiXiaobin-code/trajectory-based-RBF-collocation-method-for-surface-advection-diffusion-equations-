clc;
clear;
close;



surftype = 'Bretzel2';
load('Bretzel2_414.mat')
p0 = p;
load('Bretzel2_7216.mat')

surftype = 'Sphere';
load('Sphere_198.mat');
p0 = p;
load('Sphere_4482.mat');

surftype = 'CPD surface';
load('CPD surface_208.mat')
p0 = p;
load('CPD surface_6256.mat');


surftype = 'Torus';
load('Torus_536.mat');
p0 = p;
 load('Torus_5312.mat')



syms x y z;
switch surftype
    case 'Sphere' 
        R = 1;
        shape_f = @(x, y, z) x.^2 + y.^2 + z.^2 - R.^2;
    case 'Torus'
        R = 1;
        r = 1/3;
        shape_f = @(x, y, z) (x.^2 + y.^2 + z.^2 + R.^2 - r.^2).^2 - 4*R.^2.*(x.^2 + y.^2);
    case 'Bretzel2'
        shape_f = @(x, y, z) (-x.^4 + x.^2 - y.^2).^2 + z.^2/2 - (x.^2 + y.^2 + z.^2 + 1)/40;
    case 'CPD surface'
        shape_f = @(x, y, z) sqrt((x - 1).^2 + y.^2 + z.^2).*sqrt((x + 1).^2 + y.^2 + z.^2).*sqrt((y - 1).^2 + x.^2 + z.^2).*sqrt((y + 1).^2 + x.^2 + z.^2) - 1.1;
end


N1 = diff(shape_f, x);
N2 = diff(shape_f, y);
N3 = diff(shape_f, z);
N{1} = matlabFunction(N1./sqrt(N1^2 + N2^2 + N3^2));
N{2} = matlabFunction(N2./sqrt(N1^2 + N2^2 + N3^2));
N{3} = matlabFunction(N3./sqrt(N1^2 + N2^2 + N3^2));


switch surftype
    case 'Sphere' 
        v = @(x)  0.01*[-N{2}(x(:,1),x(:,2),x(:,3)),N{1}(x(:,1),x(:,2),x(:,3)),0*x(:,3)];
    case 'Torus'
        R = 1;
        r = 1/3;
        v = @(x) 0.01*velocityTorus(x,R,r);  
    case 'Bretzel2' 
        v = @(x)  0.01*[0*x(:,1),-N{3}(x(:,1),x(:,2),x(:,3)),N{2}(x(:,1),x(:,2),x(:,3))];
       
    case 'CPD surface'
        v = @(x)  0.01*[-N{3}(x(:,1),x(:,2),x(:,3)),0*x(:,2),N{1}(x(:,1),x(:,2),x(:,3))];
       
end
v0 = v(p0);  


% 创建图形并绘制三维曲面
fig = figure;
h = trisurf(t, p(:,1), p(:,2), p(:,3), 0*p(:,3));
hold on;
plotquiver3(p0+0.02*[N{1}(p0(:,1),p0(:,2),p0(:,3)),N{2}(p0(:,1),p0(:,2),p0(:,3)),N{3}(p0(:,1),p0(:,2),p0(:,3)),], v0);


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
% filename = sprintf('%s.png', surftype);
% print('-dpng', '-r400', filename);
frame_file = sprintf('%s/stripe_%s.png', output_folder, surftype);
exportgraphics(fig, frame_file, 'Resolution', 200);
output_folder = 'C:\Users\Admin\Desktop\fig';



