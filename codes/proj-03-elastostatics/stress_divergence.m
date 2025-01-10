function [div_sigma_x, div_sigma_y] = stress_divergence(x, y, sigma_xx, sigma_yy, sigma_xy)

% dx = x(1, 2) - x(1, 1); % x 方向步长
% dy = y(2, 1) - y(1, 1); % y 方向步长
% [sigma_xx_x, sigma_xx_y] = gradient(sigma_xx, dx, dy); % sigma_xx 的偏导数
% [sigma_xy_x, sigma_xy_y] = gradient(sigma_xy, dx, dy); % sigma_xy 的偏导数
% [sigma_yy_x, sigma_yy_y] = gradient(sigma_yy, dx, dy); % sigma_yy 的偏导数

[sigma_xx_x, sigma_xx_y] = gradient(sigma_xx, x(1, :), y(:, 1)); % 使用完整网格定义
[sigma_xy_x, sigma_xy_y] = gradient(sigma_xy, x(1, :), y(:, 1));
[sigma_yy_x, sigma_yy_y] = gradient(sigma_yy, x(1, :), y(:, 1));

% 计算散度
div_sigma_x = sigma_xx_x + sigma_xy_y; % ∂(σ_xx)/∂x + ∂(σ_xy)/∂y
div_sigma_y = sigma_xy_x + sigma_yy_y; % ∂(σ_xy)/∂x + ∂(σ_yy)/∂y
end
