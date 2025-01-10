function [div_sigma_x, div_sigma_y] = stress_divergence(x, y, Tx, R)

r = sqrt(x.^2 + y.^2);
theta = atan2(y, x);
% 若 r 在孔内，没有应力
if r < R
    sigma_xx = 0;
    sigma_yy = 0;
    sigma_xy = 0;
    return;
end

% 极坐标下应力分量
sigma_rr =  (Tx / 2) * (1 - (R^2 ./ r.^2)) + (Tx / 2) * (1 - 4 * (R^2 ./ r.^2) + 3 * (R^4 ./ r.^4)) .* cos(2 * theta);
sigma_tt =  (Tx / 2) * (1 + (R^2 ./ r.^2)) - (Tx / 2) * (1 + 3 * (R^4 ./ r.^4)) .* cos(2 * theta);
sigma_rt = -(Tx / 2) * (1 + 2 * (R^2 ./ r.^2) - 3 * (R^4 ./ r.^4)) .* sin(2 * theta);

% 转换到直角坐标系
sigma_xx = sigma_rr .* cos(theta).^2 + sigma_tt .* sin(theta).^2 - 2 * sigma_rt .* sin(theta) .* cos(theta);% 往回转theta角, 所以用-theta
sigma_yy = sigma_rr .* sin(theta).^2 + sigma_tt .* cos(theta).^2 + 2 * sigma_rt .* sin(theta) .* cos(theta);
sigma_xy = (-sigma_tt + sigma_rr) .* sin(theta) .* cos(theta) + sigma_rt .* (cos(theta).^2 - sin(theta).^2);




% % 计算散度
% [sigma_xx_x, sigma_xx_y] = gradient(sigma_xx, x(1, :), y(:, 1));
% [sigma_xy_x, sigma_xy_y] = gradient(sigma_xy, x(1, :), y(:, 1));
% [sigma_yy_x, sigma_yy_y] = gradient(sigma_yy, x(1, :), y(:, 1));
% 

% 对应力分量求梯度
dx = 0.1; % x 方向步长
dy = 0.1; % y 方向步长
[sigma_xx_x, sigma_xx_y] = gradient(sigma_xx, dx, dy); % sigma_xx 的偏导数
[sigma_xy_x, sigma_xy_y] = gradient(sigma_xy, dx, dy); % sigma_xy 的偏导数
[sigma_yy_x, sigma_yy_y] = gradient(sigma_yy, dx, dy); % sigma_yy 的偏导数


div_sigma_x = sigma_xx_x + sigma_xy_y; % ∂(σ_xx)/∂x + ∂(σ_xy)/∂y
div_sigma_y = sigma_xy_x + sigma_yy_y; % ∂(σ_xy)/∂x + ∂(σ_yy)/∂y
end
