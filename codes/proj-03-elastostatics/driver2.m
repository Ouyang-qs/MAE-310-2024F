clear; clc;
run('newone.m');

% mesh
IEN = msh.QUADS(:,1:4); % n_el * n_en
x_coor = msh.POS(:,1);
y_coor = msh.POS(:,2);

n_np = size(x_coor,1); % total number of nodal points
n_el = size(IEN,1);    % total number of elements
n_en = 4;    % number of nodes in an element (quadrilateral)
n_sd = 2;    % number of spatial dimension
niu  = 0.3;  % Possion ratio
E    = 1E6;  % Young's modulus (unit: KPa)   ??
Tx   = 10;   % traction (unit: KPa)  ??
R    = 0.5;  % radius (unit: m)
L    = 10;   % length (unit: m)
f1   = 0;    % no body force
f2   = 0;


x_coor = x_coor / L;  % 单位化
y_coor = y_coor / L;
R = R / L;

% plane stress
DD=zeros(3);
DD(1,1)=1;   DD(1,2)=niu;  
DD(2,1)=niu; DD(2,2)=1;  DD(3,3) = 0.5*(1-niu);
DD = DD * E / (1-niu^2);

% % plane strain
% lamda = niu*E/(1+niu)/(1-2*niu);
% miu   = E/2/(1+niu);
% DD = zeros(3);
% DD(1,1) = 2*miu +lamda;  DD(1,2) = lamda;  
% DD(2,1) = lamda;         DD(2,2) = 2*miu +lamda; DD(3,3) = miu;


% quadrature rule
n_int_xi  = 7;
n_int_eta = 7;
n_int     = n_int_xi * n_int_eta;
[xi, eta,  weight] = Gauss2D(n_int_xi, n_int_eta);

% Dirichlet boundary nodes
Dirichlet_BC_x = zeros(1,2);
Dirichlet_BC_y = Dirichlet_BC_x;

counter = 0;
for index = 1:size(msh.LINES,1)
    if msh.LINES(index,3)==9
        counter = counter+1;
        Dirichlet_BC_y(counter,1:2)=[msh.LINES(index,1),msh.LINES(index,2)];
    end
end

counter = 0;
for index = 1:size(msh.LINES,1)
    if msh.LINES(index,3)==10
        counter = counter+1;
        Dirichlet_BC_x(counter,1:2)=[msh.LINES(index,1),msh.LINES(index,2)];
    end
end

temp = Dirichlet_BC_x(size(Dirichlet_BC_x,1),2); % temp:temporary
Dirichlet_BC_x(size(Dirichlet_BC_x,1)+1,1) = temp;
Dirichlet_BC_x(:,2)=1; % direction1 (ii=1)

temp = Dirichlet_BC_y(size(Dirichlet_BC_y,1),2);
Dirichlet_BC_y(size(Dirichlet_BC_y,1)+1,1) = temp;
Dirichlet_BC_y(:,2)=2; % direction2 (ii=2)

Dirichlet_BC = zeros(1,2);

s1=size(Dirichlet_BC_x,1);
s2=size(Dirichlet_BC_y,1);

for index = 1:s1
    Dirichlet_BC = Dirichlet_BC_x;
end

for index = 1:s2
    Dirichlet_BC(index+s1,1) = Dirichlet_BC_y(index,1);
    Dirichlet_BC(index+s1,2) = Dirichlet_BC_y(index,2);
end

% ID array
ID = zeros(n_sd,n_np)+1;
for PP = 1 : n_np
    for index = 1 : size(Dirichlet_BC,1)
        if Dirichlet_BC(index,1) == PP
            ID( Dirichlet_BC(index,2), PP) = 0;
        end
    end
end

counter = 0;
for PP = 1 : n_np
    for ii = 1 : n_sd
        if ID(ii,PP)==0
            continue
        else 
            counter   = counter+1;
            ID(ii,PP) = counter;
        end
    end
end

n_eq = counter;

% allocate the stiffness matrix and load vector
K = spalloc(n_eq, n_eq, 16 * n_eq);  % need consideration
F = zeros(n_eq, 1);

% loop over element to assembly the matrix and vector
for ee = 1 : n_el
    x_ele = x_coor( IEN(ee, 1:n_en) );
    y_ele = y_coor( IEN(ee, 1:n_en) );

    k_ele = zeros(n_en*n_sd, n_en*n_sd); % element stiffness matrix
    f_ele = zeros(n_en*n_sd, 1);         % element load vector

    for ll = 1 : n_int
        x_l = 0.0; y_l = 0.0;
        dx_dxi = 0.0; dx_deta = 0.0;
        dy_dxi = 0.0; dy_deta = 0.0;
        for aa = 1 : n_en
            x_l = x_l + x_ele(aa) * Quad(aa, xi(ll), eta(ll));
            y_l = y_l + y_ele(aa) * Quad(aa, xi(ll), eta(ll));
            [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
            dx_dxi  = dx_dxi  + x_ele(aa) * Na_xi;
            dx_deta = dx_deta + x_ele(aa) * Na_eta;
            dy_dxi  = dy_dxi  + y_ele(aa) * Na_xi;
            dy_deta = dy_deta + y_ele(aa) * Na_eta;
        end

        detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;

        for aa = 1 : n_en
            Na = Quad(aa, xi(ll), eta(ll));
            [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
            Na_x = ( Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
            Na_y = (-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ;

            Ba = zeros(3,2); % n_sd=2;
            Ba(:,1) = [Na_x, 0, Na_y];
            Ba(:,2) = [0, Na_y, Na_x];

            for ii=1:n_sd
                pp = n_sd * (aa-1) + ii;

                if ii==1
                    f_ele(pp) = f_ele(pp) + weight(ll) * detJ * f1 * Na;
                elseif ii==2
                    f_ele(pp) = f_ele(pp) + weight(ll) * detJ * f2 * Na;
                end
                for bb = 1 : n_en
                    Nb = Quad(bb, xi(ll), eta(ll));
                    [Nb_xi, Nb_eta] = Quad_grad(bb, xi(ll), eta(ll));
                    Nb_x = ( Nb_xi * dy_deta - Nb_eta * dy_dxi) / detJ;
                    Nb_y = (-Nb_xi * dx_deta + Nb_eta * dx_dxi) / detJ;

                    Bb = zeros(3,2);
                    Bb(:,1) = [Nb_x, 0, Nb_y];
                    Bb(:,2) = [0, Nb_y, Nb_x];
                    
                    K_ab = Ba' * DD * Bb * detJ * weight(ll);

                    for jj=1:n_sd
                        qq = n_sd * (bb-1) + jj;                       
                        k_ele(pp, qq) = k_ele(pp, qq) + K_ab(ii,jj);
                    end
                end % end of bb loop
            end
        end % end of aa loop
    end % end of quadrature loop

    % assembly K,F
    for aa = 1 : n_en
        AA = IEN(ee, aa);
        for ii=1:n_sd
            PP = ID(ii, AA);
            pp = n_sd * (aa-1) + ii;
            if PP > 0
                F(PP) = F(PP) + f_ele(pp);
                for bb = 1 : n_en
                    BB = IEN(ee, bb);
                    for jj=1:n_sd
                        QQ = ID(jj, BB);
                        qq = n_sd * (bb-1) + jj;
                        if QQ > 0
                            K(PP, QQ) = K(PP, QQ) + k_ele(pp, qq);
                        else
                            % modify F with the boundary data
                            [g1,g2] = g_win( x_coor(BB), y_coor(BB) ); % g_win is Dirichlet BC
                            if jj==1
                                F(PP) = F(PP) - k_ele(pp, qq) * g1;
                            elseif jj==2
                                F(PP) = F(PP) - k_ele(pp, qq) * g2;
                            end % end of Dirichlet BC loop
                        end
                    end
                end % end of bb loop
            end
        end
    end % end of aa loop
end % end of the element loop



% Neumann boundary nodes
Neumann_BC_x = zeros(1,2);
Neumann_BC_y = Neumann_BC_x;

counter = 0;
for index = 1:size(msh.LINES,1)
    if msh.LINES(index,3)==11  % LINE with tag #11 on the right (constant x_coor)
        counter = counter+1;
        Neumann_BC_x(counter,1:2)=[msh.LINES(index,1),msh.LINES(index,2)];
    end
end

counter = 0;
for index = 1:size(msh.LINES,1)
    if msh.LINES(index,3)==8   % LINE with tag #8 on the top (constant y_coor)
        counter = counter+1;
        Neumann_BC_y(counter,1:2)=[msh.LINES(index,1),msh.LINES(index,2)];
    end
end

s1 = size(Neumann_BC_x,1);  % number of linear elements
s2 = size(Neumann_BC_y,1);

% quadrature rule
qua = 7;
[xi_1D, weight_1D] = Gauss(qua,-1,1);

% right boundary
for ss = 1:s1
    y_ele_node = zeros(2,1);
    h_ele = zeros(2,2);

    for aa = 1:2 % n_en=2 
        AA = Neumann_BC_x(ss,aa);
        y_ele_node(aa) = y_coor(AA);
    end
    % x_coor(AA) 相等

    for ll = 1 : qua
        y_l = 0.0;
        dy_dxi = 0.0;
        for aa = 1:2
            y_l    = y_l    + y_ele_node(aa) * PolyShape(1, aa, xi_1D(ll), 0);
            dy_dxi = dy_dxi + y_ele_node(aa) * PolyShape(1, aa, xi_1D(ll), 1);
        end
        
        [h1, h2] = h_win(x_coor(AA), y_l, R);

        for ii = 1 : n_sd
            for aa = 1:2
                if ii == 1
                    h_ele(ii,aa) = h_ele(ii,aa) + weight_1D(ll) * PolyShape(1, aa, xi_1D(ll), 0) * h1*Tx * dy_dxi;
                elseif ii == 2
                    h_ele(ii,aa) = h_ele(ii,aa) + weight_1D(ll) * PolyShape(1, aa, xi_1D(ll), 0) * h2*Tx * dy_dxi;
                end
            end
        end
    end
    
    for ii = 1 : n_sd
        for aa = 1:2
            AA = Neumann_BC_x(ss,aa);
            PP = ID(ii, AA);
            if (PP > 0)   % 排除同时处在两个边界上的点 ( g and h )
                F(PP) = F(PP) + h_ele(ii,aa);
            end
        end
    end
end

% top boundary
for ss = 1:s2
    x_ele_node = zeros(2,1);
    h_ele = zeros(2,2);

    for aa = 1:2 % n_en=2 
        AA = Neumann_BC_y(ss,aa);
        x_ele_node(aa) = x_coor(AA);
    end
    % y_coor(AA) 相等

    for ll = 1 : qua
        x_l = 0.0;
        dx_dxi = 0.0;
        for aa = 1:2
            x_l    = x_l    + x_ele_node(aa) * PolyShape(1, aa, xi_1D(ll), 0);
            dx_dxi = dx_dxi + x_ele_node(aa) * PolyShape(1, aa, xi_1D(ll), 1);
        end

        [h1, h2] = h_win(x_l, y_coor(AA), R);

        for ii = 1 : n_sd
            for aa = 1:2
                if ii == 1
                    h_ele(ii,aa) = h_ele(ii,aa) + weight_1D(ll) * PolyShape(1, aa, xi_1D(ll), 0) * h1*Tx * dx_dxi;
                elseif ii == 2
                    h_ele(ii,aa) = h_ele(ii,aa) + weight_1D(ll) * PolyShape(1, aa, xi_1D(ll), 0) * h2*Tx * dx_dxi;
                end
            end
        end
    end

    for ii = 1 : n_sd
        for aa = 1:2
            AA = Neumann_BC_y(ss,aa);
            PP = ID(ii, AA);
            if (PP > 0)   % 排除同时处在两个边界上的点 ( g and h )
                F(PP) = F(PP) + h_ele(ii,aa);
            end
        end
    end
end

% solve the stiffness matrix
dn = K \ F;

% insert dn back into the vector for all nodes
disp = zeros(n_np, n_sd);

for nn = 1 : n_np
    for ii = 1 : n_sd
        index = ID(ii,nn);
        if index > 0
            disp(nn,ii) = dn(index);
        else
            % modify disp with the g data.
            [g1,g2] = g_win( x_coor(nn), y_coor(nn) ); % g_win is Dirichlet BC
            if ii==1
                disp(nn,ii) = g1;
            elseif ii==2
                disp(nn,ii) = g2;
            end         
        end
    end
end

disp = disp * L;
x_coor = x_coor * L;
y_coor = y_coor * L;

% save the solution vector and number of elements to disp with name
% ELASTO.mat
save("ELASTO2", "disp","x_coor","y_coor","");



% visualization
disp_x = disp(:,1);
disp_y = disp(:,2);

IEN_tri = zeros(1,1);
for ee = 1:size(IEN,1)
    IEN_tri(ee*2-1,1) = IEN(ee,1);
    IEN_tri(ee*2-1,2) = IEN(ee,2);
    IEN_tri(ee*2-1,3) = IEN(ee,3);
    IEN_tri(ee*2,1) = IEN(ee,1);
    IEN_tri(ee*2,2) = IEN(ee,3);
    IEN_tri(ee*2,3) = IEN(ee,4);
end

% Displacement
figure(1)
% hold on;
trisurf(IEN_tri, x_coor, y_coor, disp_x);

colormap jet
shading interp
colorbar;
% axis equal;

title('Displacement Component x');
xlabel('x');
ylabel('y');


figure(2)
% hold on;
trisurf(IEN_tri, x_coor, y_coor, disp_y);

colormap jet
shading interp
colorbar;
% axis equal;

title('Displacement Component y');
xlabel('x');
ylabel('y');




% Strain and Stress
n = n_np;

[ux_dx, uy_dy, ux_dy, uy_dx] = derivative(x_coor, y_coor, disp, IEN);

epsilon_vect = zeros(n,3);
for i = 1:n
    epsilon_vect(i,1) = ux_dx(i);
    epsilon_vect(i,2) = uy_dy(i);
    epsilon_vect(i,3) = uy_dx(i) + ux_dy(i);
end

sigma_vect=zeros(n,3);
for i=1:n
    sigma_vect(i,1:3) = DD * epsilon_vect(i,1:3)';
end

epsilon11 = epsilon_vect(:,1);
epsilon22 = epsilon_vect(:,2);
epsilon12 = epsilon_vect(:,3);

sigma11 = sigma_vect(:,1);
sigma22 = sigma_vect(:,2);
sigma12 = sigma_vect(:,3);



figure(3)
% hold on;
trisurf(IEN_tri, x_coor, y_coor, epsilon11);
colormap jet
shading interp
colorbar;
% axis equal;

title('Strain Component x');
xlabel('x');
ylabel('y');


figure(4)
% hold on;
trisurf(IEN_tri, x_coor, y_coor, sigma11);

colormap jet
shading interp
colorbar;
% axis equal;

title('Stress Component x');
xlabel('x');
ylabel('y');

