clear all; clc;

n_sd = 2;   % number of spatial dimension
niu = 0.3;  % Possion ratio
E   = 100;  % Young's modulus (unit:GPa)

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

% manufacture solution
% exact_x = @(x,y) x*(1-x)*y*(1-y); 
% exact_y = @(x,y) x*(1-x)*y*(1-y);

exact_x = @(x,y) x^2+x*y;
exact_y = @(x,y) x+2*y;


% f1 = @(x,y) E/(1-niu^2) * (2*y*(1-y)-niu*(1-2*x)*(1-2*y)) + E/(2+2*niu)*(2*x*(1-x)-(1-2*x)*(1-2*y));
% f2 = @(x,y) E/(1-niu^2) * (2*x*(1-x)-niu*(1-2*x)*(1-2*y)) + E/(2+2*niu)*(2*y*(1-y)-(1-2*x)*(1-2*y));

f1 = @(x,y) E/(1-niu^2) * (-2);
f2 = @(x,y) E/(1-niu^2) * (1+niu) / (-2);


% quadrature rule
n_int_xi  = 5;
n_int_eta = 5;
n_int     = n_int_xi * n_int_eta;
[xi, eta,  weight] = Gauss2D(n_int_xi, n_int_eta);
[xi_1D, weight_1D] = Gauss(n_int,0,1);


% mesh generation
n_en   = 4;               % number of nodes in an element (quadrilateral)
n_el_x = 30;              % number of elements in x-dir
n_el_y = 30;              % number of elements in y-dir
n_el   = n_el_x * n_el_y; % total number of elements

n_np_x = n_el_x + 1;      % number of nodal points in x-dir
n_np_y = n_el_y + 1;      % number of nodal points in y-dir
n_np   = n_np_x * n_np_y; % total number of nodal points

x_coor = zeros(n_np, 1);
y_coor = x_coor;

hx = 1.0 / n_el_x;        % mesh size in x-dir
hy = 1.0 / n_el_y;        % mesh size in y-dir

% generate the nodal coordinates
for ny = 1 : n_np_y
    for nx = 1 : n_np_x
        index = (ny-1)*n_np_x + nx; % nodal index
        x_coor(index) = (nx-1) * hx;
        y_coor(index) = (ny-1) * hy;
    end
end

% IEN array
IEN = zeros(n_el, n_en);
for ex = 1 : n_el_x
    for ey = 1 : n_el_y
        ee = (ey-1) * n_el_x + ex; % element index
        IEN(ee, 1) = (ey-1) * n_np_x + ex;
        IEN(ee, 2) = (ey-1) * n_np_x + ex + 1;
        IEN(ee, 3) =  ey    * n_np_x + ex + 1;
        IEN(ee, 4) =  ey    * n_np_x + ex;
    end
end

% ID array
ID = zeros(n_sd,n_np);
counter = 0;
for ny = 2 : n_np_y - 1
    for nx = 2 : n_np_x - 1
        index = (ny-1)*n_np_x + nx;
        counter = counter + 1;
        ID(index*2-1) = counter;
        counter = counter + 1;
        ID(index*2) = counter;
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
                    f_ele(pp) = f_ele(pp) + weight(ll) * detJ * f1(x_l, y_l) * Na;
                elseif ii==2
                    f_ele(pp) = f_ele(pp) + weight(ll) * detJ * f2(x_l, y_l) * Na;
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

    for aa = 1 : n_en
        index_1 = IEN(ee, aa);
        for ii=1:n_sd
            PP = ID(ii, index_1);
            pp = n_sd * (aa-1) + ii;
            if PP > 0
                F(PP) = F(PP) + f_ele(pp);
                for bb = 1 : n_en
                    index_2 = IEN(ee, bb);
                    for jj=1:n_sd
                        QQ = ID(jj, index_2);
                        qq = n_sd * (bb-1) + jj;
                        if QQ > 0
                            K(PP, QQ) = K(PP, QQ) + k_ele(pp, qq);
                        else
                            % modify F with the boundary data
                            [g1,g2] = g( x_coor(index_2), y_coor(index_2) ); % g is a function
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
            [g1,g2] = g( x_coor(nn), y_coor(nn) ); % g is a function
            if ii==1
                disp(nn,ii) = g1;
            elseif ii==2
                disp(nn,ii) = g2;
            end         
        end
    end
end

% save the solution vector and number of elements to disp with name
% ELASTO.mat
save("ELASTO", "disp", "n_el_x", "n_el_y");


% EOF