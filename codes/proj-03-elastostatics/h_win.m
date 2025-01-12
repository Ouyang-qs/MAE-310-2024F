function [h1,h2] = h_win(x, y, R)

[sigma_xx, sigma_yy, sigma_xy] = stress_cartesian(x, y, R);

[n1, n2] = unit_normal_verctor(x,y);

h1 = sigma_xx*n1 + sigma_xy*n2;
h2 = sigma_xy*n1 + sigma_yy*n2;

end

function [n1,n2] = unit_normal_verctor(x,y)
% Neumann boundary conditions
% nothing to do with the arc (no hi on it)

if x==0 && y~=0
    n1 = -1;
    n2 = 0;
elseif y==0 && x~=0
    n1 = 0;
    n2 = -1;
elseif y==1 && x==1
    n1 = 1;
    n2 = 1;
elseif x==1 && y~=1
    n1 = 1;
    n2 = 0;
elseif x~=1 && y==1
    n1 = 0;
    n2 = 1;
elseif x==0 && y==0
    n1 = 0;
    n2 = 0;
end

end


% EOF
