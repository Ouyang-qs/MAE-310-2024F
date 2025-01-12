function [dux_dx, duy_dy, dux_dy, duy_dx] = derivative(x_coor, y_coor, disp, IEN)
% x_coor, y_coor: 节点的坐标
% disp位移: 第一列为x方向位移，第二列为y方向位移
% IEN: 元素节点连接信息

n_np = length(x_coor);      % 总节点数
n_el = size(IEN, 1);        % 元素数量

dux_dx = zeros(n_np, 1);    % 存储 u_x 对 x 的偏导数
duy_dy = dux_dx;            % 存储 u_y 对 y 的偏导数
dux_dy = dux_dx;            % 存储 u_x 对 y 的偏导数
duy_dx = dux_dx;            % 存储 u_y 对 x 的偏导数


% 计算每个节点的导数
for i = 1:n_np
    % 通过 IEN 数组获取邻接节点，遍历每个元素，找到相邻的节点
    dux_dx_sum = 0;
    duy_dy_sum = 0;
    dux_dy_sum = 0;
    duy_dx_sum = 0;

    neighbor_count_x = 0;
    neighbor_count_y = 0;

    for ee = 1:n_el
        % 检查当前节点是否在元素内
        if any(IEN(ee, :) == i)  % 如果当前节点i属于元素ee
            % 获取该元素内的其他节点
            neighbors = setdiff(IEN(ee, :), i);  % 获取当前节点的邻接节点
            for AA = neighbors
                % 计算相邻节点的差分
                if (x_coor(AA) - x_coor(i)) > 1E-4
                    dux_dx_sum = dux_dx_sum + (disp(AA, 1) - disp(i, 1)) / (x_coor(AA) - x_coor(i));

                    duy_dx_sum = duy_dx_sum + (disp(AA, 2) - disp(i, 2)) / (x_coor(AA) - x_coor(i));

                    neighbor_count_x = neighbor_count_x + 1;
                end

                if (y_coor(AA) - y_coor(i)) > 1E-5
                    duy_dy_sum = duy_dy_sum + (disp(AA, 2) - disp(i, 2)) / (y_coor(AA) - y_coor(i));

                    dux_dy_sum = dux_dy_sum + (disp(AA, 1) - disp(i, 1)) / (y_coor(AA) - y_coor(i));

                    neighbor_count_y = neighbor_count_y + 1;
                end
            end
        end
    end

    if neighbor_count_x > 0
        dux_dx(i) = dux_dx_sum / neighbor_count_x;  % 平均值

        duy_dx(i) = duy_dx_sum / neighbor_count_x;
    end

    if neighbor_count_y > 0
        duy_dy(i) = duy_dy_sum / neighbor_count_y;

        dux_dy(i) = dux_dy_sum / neighbor_count_y;
    end
end

end

