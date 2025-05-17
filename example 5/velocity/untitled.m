% 定义两个漩涡中心的位置
x_c1 = -0.6;
y_c1 = 0;
x_c2 = 0.6;
y_c2 = 0;

% 定义网格
[x, y] = meshgrid(linspace(-2, 2, 40), linspace(-2, 2, 40));
z = zeros(size(x));  % 在 xy 平面上绘制

% 计算第一个漩涡中心的速度向量
vx1 = y - y_c1;
vy1 = -(x - x_c1);
vz1 = zeros(size(x));

% 计算第二个漩涡中心的速度向量
vx2 = y - y_c2;
vy2 = -(x - x_c2);
vz2 = zeros(size(x));

% 将两个速度向量叠加
vx = vx1 + vx2;
vy = vy1 + vy2;
vz = vz1 + vz2;

% 绘制速度向量场
figure;
quiver3(x, y, z, vx, vy, vz, 'AutoScaleFactor', 1.5);
axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Velocity Field for Rotation around Two Points');
grid on;
