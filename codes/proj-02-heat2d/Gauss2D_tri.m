function [xi, eta, w] = Gauss2D_tri(N1, N2)

% preallocation
xi = zeros(N1*N2,1);
eta = xi;
w = xi;

% generate 1D rule
[x1, w1] = Gauss(N1, -1, 1);  % [-1, 1] 高斯节点和权重
[x2, w2] = Gauss(N2, -1, 1);  % [-1, 1] 高斯节点和权重

for ii = 1 : N1
    for jj = 1 : N2
        % 对每个高斯点进行坐标变换
        index = (jj-1)*N1 + ii;
        
        % 高斯节点变换
        x = (x1(ii) + 1) / 2;  % 从 [-1, 1] 到 [0, 1]
        y = (x2(jj) + 1) / 2 * (1 - x);  % 从 [-1, 1] 到 [0, 1-x]

        xi(index) = x;  % 记录新的 x 坐标
        eta(index) = y;  % 记录新的 y 坐标
        
        % 权重变换
        w(index) = w1(ii) * w2(jj) * (1 - x)/4;  % 根据雅可比变换计算权重
    end
end

% EOF
