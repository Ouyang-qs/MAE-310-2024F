function [g1,g2] = g(x,y)

% exact_x = @(x,y) x^2+x*y;
% exact_y = @(x,y) x+2*y;
% assume all edges are controled by Dirchilet boundary conditions.
index=0;
if x==0
    g1 = 0;
    g2 = 2*y;
else
    index = index +1;
end

if y==0
    g1 = x^2;
    g2 = x;
else
    index = index +1;
end

if x==1
    g1 = 1+y;
    g2 = 1+2*y;
else
    index = index +1;
end

if y==1
    g1 = x^2+x;
    g2 = x+2;
else
    index = index +1;
end

if index==4
    error('Error: the displacement of this point is not determined by Dirichlet BC.');
end




% EOF