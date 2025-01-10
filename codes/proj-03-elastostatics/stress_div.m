function [div_sigma_x_value, div_sigma_y_value] = stress_div(x_value, y_value, T, R) 

n = 100; % 网格大小
x_range = linspace(0, 1, n);           % x 方向范围
y_range = linspace(0, 1, n);           % y 方向范围
[x, y] = meshgrid(x_range, y_range);   % 生成网格

% 极坐标 r,theta
r = sqrt(x.^2 + y.^2); 
theta = atan2(y, x);   

if r < R
    div_sigma_x_value = 0;
    div_sigma_y_value = 0;
    return;
end

% 极坐标下的应力分量
sigma_rr = T / 2 * (1 - R^2 ./ r.^2) + T / 2 * (1 - 4 * R^2 ./ r.^2 + 3 * R^4 ./ r.^4) .* cos(2 * theta);
sigma_tt = T / 2 * (1 + R^2 ./ r.^2) - T / 2 * (1 + 3 * R^4 ./ r.^4) .* cos(2 * theta);
sigma_rt = -T / 2 * (1 + 2 * R^2 ./ r.^2 - 3 * R^4 ./ r.^4) .* sin(2 * theta);

% 转换到直角坐标系
sigma_xx =  sigma_rr .* cos(theta).^2 + sigma_tt .* sin(theta).^2 - 2 * sigma_rt .* sin(theta) .* cos(theta);
sigma_yy =  sigma_rr .* sin(theta).^2 + sigma_tt .* cos(theta).^2 + 2 * sigma_rt .* sin(theta) .* cos(theta);
sigma_xy = (sigma_rr - sigma_tt) .* sin(theta) .* cos(theta) + sigma_rt .* (cos(theta).^2 - sin(theta).^2);

% 应力分量的散度
dx = x(1, 2) - x(1, 1); % x 方向步长
dy = y(2, 1) - y(1, 1); % y 方向步长

[sigma_xx_x, sigma_xx_y] = gradient(sigma_xx, dx, dy); % sigma_xx 的偏导数
[sigma_xy_x, sigma_xy_y] = gradient(sigma_xy, dx, dy); % sigma_xy 的偏导数
[sigma_yy_x, sigma_yy_y] = gradient(sigma_yy, dx, dy); % sigma_yy 的偏导数

div_sigma_x = sigma_xx_x + sigma_xy_y; % ∂(σ_xx)/∂x + ∂(σ_xy)/∂y
div_sigma_y = sigma_xy_x + sigma_yy_y; % ∂(σ_xy)/∂x + ∂(σ_yy)/∂y


% 带入具体 (x, y) 坐标值
% 在网格中找到最近的索引
[~, ix] = min(abs(x(1, :) - x_value));
[~, iy] = min(abs(y(:, 1) - y_value));

% 获取对应散度值
div_sigma_x_value = div_sigma_x(iy, ix); % div(sigma_x) 在指定点的值
div_sigma_y_value = div_sigma_y(iy, ix); % div(sigma_y) 在指定点的值

end
