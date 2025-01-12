x = [1/30, 1/50, 1/70, 1/90];  % h

y=[2.028602e-04, 7.302967e-05,3.726003e-05, 2.254002e-05]; % L2
   
z=[0.019245, 0.0115470, 0.008247, 0.006415]; % H1


% 绘制 log-log 图
figure;
loglog(x, y, 'o', 'MarkerSize', 8, 'LineWidth', 2);

hold on;
loglog(x, z, '*', 'MarkerSize', 8, 'LineWidth', 2);

xlabel('log(h)');
ylabel('log(error)');
grid on;
title('Log-Log Plot');


% 计算 log-log 图的斜率
% 对 x 和 y 进行对数转换
log_x = log10(x);
log_y = log10(y);
log_z = log10(z);

% polyfit 拟合线性关系，1 表示直线拟合
[p, S] = polyfit(log_x, log_y, 1);  % p(1) 是斜率，p(2) 是截距
[q, Q] = polyfit(log_x, log_z, 1);


% 输出拟合的斜率
disp(['L2斜率是：', num2str(p(1))]);
disp(['H1斜率是：', num2str(q(1))]);

% 绘制拟合直线
hold on;
loglog(x, 10.^(polyval(p, log_x)), 'r-', 'LineWidth', 2);
legend('Data points', 'Fitted line');
hold off;

hold on;
loglog(x, 10.^(polyval(q, log_x)), 'r-', 'LineWidth', 2);
legend('Data points', 'Fitted line');
hold off;
