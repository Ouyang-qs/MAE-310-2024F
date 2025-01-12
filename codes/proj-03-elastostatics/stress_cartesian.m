function [sigma_xx, sigma_yy, sigma_xy] = stress_cartesian(x, y, R)


% (x,y):直角坐标系下点的坐标 （0，1）

% Tx: 远场拉伸应力 = 1 （单位化）
% R:  圆孔的半径 = R / L

Tx = 1;
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
% 往回转theta角, 所以用-theta
sigma_xx =  0.5*(sigma_rr + sigma_tt) + 0.5*(sigma_rr - sigma_tt) * cos(2*theta) - sigma_rt * sin(2*theta);
sigma_yy =  0.5*(sigma_rr + sigma_tt) - 0.5*(sigma_rr - sigma_tt) * cos(2*theta) + sigma_rt * sin(2*theta);
sigma_xy = -0.5*(sigma_rr - sigma_tt) * sin(2*theta)  + sigma_rt  * cos(2*theta);

end
