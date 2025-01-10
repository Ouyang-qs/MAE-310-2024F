function [sigma_xx, sigma_yy, sigma_xy] = stress_cartesian(x, y, Tx, R)
% (x,y):直角坐标系下点的坐标
% Tx: 远场拉伸应力
% R:  圆孔的半径

% 极坐标 r 和 theta
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

end
