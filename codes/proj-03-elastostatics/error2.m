% 输入坐标点
x = [2/5, 2/10, 2/20];  % h

y = [5.958810877873276e+03, 6.210039508102203e+03, 6.299067008419257e+03, 6.333e+03]; % L2

% 绘制 log-log 图
figure;
loglog(x, y, 'o', 'MarkerSize', 8, 'LineWidth', 2);

xlabel('log(h)');
ylabel('log(error)');
grid on;
title('Log-Log Plot');


% 计算 log-log 图的斜率
% 对 x 和 y 进行对数转换
log_x = log10(x);
log_y = log10(y);

% polyfit 拟合线性关系，1 表示直线拟合
[p, S] = polyfit(log_x, log_y, 1);  % p(1) 是斜率，p(2) 是截距


% 输出拟合的斜率
disp(['L2斜率是：', num2str(p(1))]);

% 绘制拟合直线
hold on;
loglog(x, 10.^(polyval(p, log_x)), 'r-', 'LineWidth', 2);
legend('Data points', 'Fitted line');
hold off;

