function [g1,g2] = g_win(x,y)
% Dirchilet boundary conditions for window problem

if x==0
    g1 = 0;
    g2 = 0;
end

if y==0
    g1 = 0;
    g2 = 0;
end

end
