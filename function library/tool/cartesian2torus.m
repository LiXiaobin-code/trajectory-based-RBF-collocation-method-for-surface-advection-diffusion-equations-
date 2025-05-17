function [theta, phi] = cartesian_to_torus(X, Y, Z, R, r)
% CARTESIAN_TO_TORUS 将笛卡尔坐标转换为环面参数（修正版本）
% 输入：
%   X, Y, Z : 三维笛卡尔坐标
%   R       : 环面主半径
%   r       : 管道半径
% 输出：
%   theta   : 绕管道方向的角参数 [0, 2π]
%   phi     : 绕环面中心轴的角参数 [0, 2π]

% 步骤1：计算phi（绕中心轴的角度）
phi = atan2(Y, X);

% 步骤2：计算投影到中心平面后的径向偏移
rho = sqrt(X.^2 + Y.^2) - R;

% 步骤3：归一化局部坐标系（关键修正：引入r）
rho_normalized = rho / r;
z_normalized = Z / r;

% 步骤4：计算theta（绕管道方向的角度）
theta = atan2(z_normalized, rho_normalized);

% 角度规范化到[0, 2π)
theta = mod(theta, 2*pi);
phi = mod(phi, 2*pi);
end