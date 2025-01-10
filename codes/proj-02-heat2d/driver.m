clear all; clc; close all

kappa = 1.0; % conductivity

% exact solution
exact = @(x,y) x*(1-x)*y*(1-y);
exact_x = @(x,y) (1-2*x)*y*(1-y);
exact_y = @(x,y) x*(1-x)*(1-2*y);

f = @(x,y) 2.0*kappa*x*(1-x) + 2.0*kappa*y*(1-y); % source term

% quadrature rule
n_int_xi  = 3;
n_int_eta = 3;
n_int     = n_int_xi * n_int_eta;
[xi, eta, weight] = Gauss2D(n_int_xi, n_int_eta);

L2 = zeros(12,1);
H1 = L2;
log_h = L2;

for iii=1:12
    % mesh generation
    n_en   = 4;               % number of nodes in an element
    n_el_x = 10*iii;              % number of elements in x-dir
    n_el_y = 10*iii;              % number of elements in y-dir
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
    ID = zeros(n_np,1);
    counter = 0;
    for ny = 2 : n_np_y - 1
        for nx = 2 : n_np_x - 1
            index = (ny-1)*n_np_x + nx;
            counter = counter + 1;
            ID(index) = counter;
        end
    end

    n_eq = counter;

    LM = ID(IEN);

    % allocate the stiffness matrix and load vector
    K = spalloc(n_eq, n_eq, 9 * n_eq);
    F = zeros(n_eq, 1);

    % loop over element to assembly the matrix and vector
    for ee = 1 : n_el
        x_ele = x_coor( IEN(ee, 1:n_en) );
        y_ele = y_coor( IEN(ee, 1:n_en) );

        k_ele = zeros(n_en, n_en); % element stiffness matrix
        f_ele = zeros(n_en, 1);    % element load vector

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
                Na_x = (Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
                Na_y = (-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ;

                f_ele(aa) = f_ele(aa) + weight(ll) * detJ * f(x_l, y_l) * Na;

                for bb = 1 : n_en
                    Nb = Quad(bb, xi(ll), eta(ll));
                    [Nb_xi, Nb_eta] = Quad_grad(bb, xi(ll), eta(ll));
                    Nb_x = (Nb_xi * dy_deta - Nb_eta * dy_dxi) / detJ;
                    Nb_y = (-Nb_xi * dx_deta + Nb_eta * dx_dxi) / detJ;

                    k_ele(aa, bb) = k_ele(aa,bb) + weight(ll) * detJ * kappa * (Na_x * Nb_x + Na_y * Nb_y);
                end % end of bb loop
            end % end of aa loop
        end % end of quadrature loop

        for aa = 1 : n_en
            PP = LM(ee, aa);
            if PP > 0
                F(PP) = F(PP) + f_ele(aa);

                for bb = 1 : n_en
                    QQ = LM(ee, bb);
                    if QQ > 0
                        K(PP, QQ) = K(PP, QQ) + k_ele(aa, bb);
                    else
                        % modify F with the boundary data
                        % here we do nothing because the boundary data g is zero or
                        % homogeneous
                    end
                end
            end
        end
    end

    % solve the stiffness matrix
    dn = K \ F;

    % insert dn back into the vector for all nodes
    disp = zeros(n_np, 1);

    for ii = 1 : n_np
        index = ID(ii);
        if index > 0
            disp(ii) = dn(index);
        else
            % modify disp with the g data. Here it does nothing because g is zero
        end
    end

    % save the solution vector and number of elements to disp with name
    % HEAT.mat
    save("HEAT", "disp", "n_el_x", "n_el_y");

    % calculate the error

    L2_top = 0.0; L2_bot = 0.0; H1_top = 0.0; H1_bot = 0.0;

    for ee = 1 : n_el
        x_ele = x_coor( IEN(ee, :) );
        y_ele = y_coor( IEN(ee, :) );
        u_ele = disp( IEN(ee, :) );

        for ll = 1 : n_int
            x_l = 0.0; y_l = 0.0;
            dx_dxi = 0.0; dx_deta = 0.0;
            dy_dxi = 0.0; dy_deta = 0.0;
            uh = 0.0; uh_xi = 0.0; uh_eta = 0.0;

            for aa = 1 : n_en
                x_l = x_l + x_ele(aa) * Quad(aa, xi(ll), eta(ll));
                y_l = y_l + y_ele(aa) * Quad(aa, xi(ll), eta(ll));
                [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
                dx_dxi  = dx_dxi  + x_ele(aa) * Na_xi;
                dx_deta = dx_deta + x_ele(aa) * Na_eta;
                dy_dxi  = dy_dxi  + y_ele(aa) * Na_xi;
                dy_deta = dy_deta + y_ele(aa) * Na_eta;

                uh      = uh      + u_ele(aa) * Quad(aa, xi(ll), eta(ll));
                uh_xi   = uh_xi   + u_ele(aa) * Na_xi;
                uh_eta  = uh_eta  + u_ele(aa) * Na_eta;

            end

            detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;

            L2_top = L2_top + weight(ll) * detJ * (uh - exact(x_l, y_l))^2;
            L2_bot = L2_bot + weight(ll) * detJ * exact(x_l, y_l)^2;

            uh_x = (uh_xi * dy_deta - uh_eta * dy_dxi) / detJ;
            uh_y = (-uh_xi * dx_deta + uh_eta * dx_dxi) / detJ;
            H1_top = H1_top + weight(ll) * detJ * ( ( uh_x - exact_x(x_l, y_l) )^2 + ( uh_y - exact_y(x_l, y_l) )^2 );
            H1_bot = H1_bot + weight(ll) * detJ * ( exact_x(x_l, y_l)^2 + exact_y(x_l, y_l)^2 );
        end
    end

    L2_top = sqrt(L2_top); L2_bot = sqrt(L2_bot);
    H1_top = sqrt(H1_top); H1_bot = sqrt(H1_bot);

    L2_error = L2_top / L2_bot;
    H1_error = H1_top / H1_bot;
    
    L2(iii) = L2_error;
    H1(iii) = H1_error;
    log_h(iii)  = log(sqrt(hx*hy));
end

log_L2 = log(L2);
log_H1 = log(H1);

plot(log_h,log_L2,'b','LineWidth',2);
hold on
plot(log_h,log_H1,'r','LineWidth',2);

title('Quadrilateral element')
xlabel('log(h)')
ylabel('log(error)')
legend('L2 Error', 'H1 Error', 'Location', 'best')
grid on


% EOF