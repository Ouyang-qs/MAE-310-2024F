clear all; clc; clf;

load("ELASTO.mat");

hh_x = 1.0 / n_el_x;
hh_y = 1.0 / n_el_y;

n_np_x = n_el_x + 1;
n_np_y = n_el_y + 1;

disp_x=disp(:,1);
disp_y=disp(:,2);

figure(2)
[X, Y] = meshgrid(0 : hh_x : 1, 0 : hh_y : 1);
Z = reshape(disp_x, n_np_x, n_np_y)';
surf(X, Y, Z);

shading interp

az = -61;
el = 20;
view(az,el);
title('displacement in x')
xlabel('x')
ylabel('y')
colormap jet;
colorbar;
%
% figure(2)
% [X, Y] = meshgrid(0 : hh_x : 1, 0 : hh_y : 1);
% Z = reshape(disp_y, n_np_x, n_np_y)';
% surf(X, Y, Z);
% 
% shading interp
% 
% az = -61;
% el = 20;
% view(az,el);
% 
% title('displacement in y')
% xlabel('x')
% ylabel('y')

% EOF