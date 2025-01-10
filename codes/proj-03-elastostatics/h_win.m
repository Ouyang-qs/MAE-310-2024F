function [h1,h2] = h_win(x, y, Tx, R)

[sigma_xx, sigma_yy, sigma_xy] = stress_cartesian(x, y, Tx, R);

[n1, n2] = unit_normal_verctor(x,y);

h1 = sigma_xx*n1 + sigma_xy*n2;
h2 = sigma_xy*n1 + sigma_yy*n2;

end

function [n1,n2] = unit_normal_verctor(x,y)
% Neumann boundary conditions
% nothing to do with the arc (no hi on it)
index=0;
if x==-1
    n1 = -1;
    n2 = 0;
else
    index = index +1;
end

if y==-1
    n1 = 0;
    n2 = -1;
else
    index = index +1;
end

if x==1
    n1 = 1;
    n2 = 0;
else
    index = index +1;
end

if y==1
    n1 = 0;
    n2 = 1;
else
    index = index +1;
end

if index==4
    error('Error: the displacement of this point is not determined by Neumann BC.');
end

end


% EOF
