clear all; clc; clf;  % clean the memory, screen, and figure

% Problem definition
f = @(x) -20*x.^3; % f(x) is the source
g = 1.0;           % u    = g  at x = 1
h = 0.0;           % -u,x = h  at x = 0

% Setup the mesh
pp   = 3;              % polynomial degree
n_en = pp + 1;         % number of element or local nodes

El2 = zeros(8,1);
Eh1 = El2;
hh  = zeros(8,1);
for s = 1:8                % n_el = 2,4...16
    n_el = 2*s;            % number of elements
    n_np = n_el * pp + 1;  % number of nodal points
    n_eq = n_np - 1;       % number of equations
    n_int = 5;

    hh(s) = 1.0 / (n_np - 1); % space between two adjacent nodes
    x_coor = 0 : hh(s) : 1;   % nodal coordinates for equally spaced nodes

    IEN = zeros(n_el, n_en);

    for ee = 1 : n_el
        for aa = 1 : n_en
            IEN(ee, aa) = (ee - 1) * pp + aa;
        end
    end

    % Setup the ID array for the problem
    ID = 1 : n_np;
    ID(end) = 0;

    % Setup the quadrature rule
    [xi, weight] = Gauss(n_int, -1, 1);

    % allocate the stiffness matrix
    K = spalloc(n_eq, n_eq, (2*pp+1)*n_eq);
    F = zeros(n_eq, 1);

    % Assembly of the stiffness matrix and load vector
    for ee = 1 : n_el
        k_ele = zeros(n_en, n_en); % allocate a zero element stiffness matrix
        f_ele = zeros(n_en, 1);    % allocate a zero element load vector

        x_ele = x_coor(IEN(ee,:)); % x_ele(aa) = x_coor(A) with A = IEN(aa, ee)

        % quadrature loop
        for qua = 1 : n_int
            dx_dxi = 0.0;
            x_l = 0.0;
            for aa = 1 : n_en
                x_l    = x_l    + x_ele(aa) * PolyShape(pp, aa, xi(qua), 0);
                dx_dxi = dx_dxi + x_ele(aa) * PolyShape(pp, aa, xi(qua), 1);
            end
            dxi_dx = 1.0 / dx_dxi;

            for aa = 1 : n_en
                f_ele(aa) = f_ele(aa) + weight(qua) * PolyShape(pp, aa, xi(qua), 0) * f(x_l) * dx_dxi;
                for bb = 1 : n_en
                    k_ele(aa, bb) = k_ele(aa, bb) + weight(qua) * PolyShape(pp, aa, xi(qua), 1) * PolyShape(pp, bb, xi(qua), 1) * dxi_dx;
                end
            end
        end

        % Assembly of the matrix and vector based on the ID or LM data
        for aa = 1 : n_en
            P = ID(IEN(ee,aa));
            if(P > 0)
                F(P) = F(P) + f_ele(aa);
                for bb = 1 : n_en
                    Q = ID(IEN(ee,bb));
                    if(Q > 0)
                        K(P, Q) = K(P, Q) + k_ele(aa, bb);
                    else
                        F(P) = F(P) - k_ele(aa, bb) * g; % handles the Dirichlet boundary data
                    end
                end
            end
        end
    end

    % ee = 1 F = NA(0)xh
    F(ID(IEN(1,1))) = F(ID(IEN(1,1))) + h;

    % Solve Kd = F equation
    d_temp = K \ F;
    
    % restart = 10000;
    % maxit   = 10000;
    % tol     = 1E-2;
    % d2 = gmres(K,F,restart,tol,maxit);
    disp  = [d_temp; g];
    % disp2 = [d2; g];
    % diff = disp-disp2;

    % plot(x_coor,diff,'b','LineWidth',1);

    
    %% Postprocessing: visualization

    %plot(x_coor, disp, '--r','LineWidth',3);

    %x_sam = 0 : 0.01 : 1;
    %y_sam = x_sam.^5;
    %hold on;
    %plot(x_sam, y_sam, '-k', 'LineWidth', 3);

    n_sam = 20;
    xi_sam = -1 : (2/n_sam) : 1;
    % xi_sam = linspace(-1,1,n_sam+1);

    x_sam = zeros(n_el * n_sam + 1, 1);
    y_sam = x_sam; % store the exact solution value at sampling points
    u_sam = x_sam; % store the numerical solution value at sampling pts
    y_x_sam = x_sam;
    u_x_sam = x_sam;

    el2_n = zeros(n_el,1); % numerator
    el2_d = el2_n;         % denominator
    El2_n = 0;
    El2_d = 0;

    eh1_n = zeros(n_el,1); % numerator
    eh1_d = eh1_n;         % denominator
    Eh1_n = 0;
    Eh1_d = 0;

    for ee = 1 : n_el
        x_ele = x_coor( IEN(ee, :) );
        u_ele = disp( IEN(ee, :) );

        if ee == n_el
            n_sam_end = n_sam+1;
        else
            n_sam_end = n_sam;
        end

        %quadrature loop
        for qua = 1 : n_int
            dx_dxi = 0.0;
            x_l    = 0.0;
            u_l    = 0.0;
            u_l_x  = 0.0;



            for aa = 1 : n_en
                x_l    = x_l    + x_ele(aa) * PolyShape(pp, aa, xi(qua), 0); % x for quadrature rule
                dx_dxi = dx_dxi + x_ele(aa) * PolyShape(pp, aa, xi(qua), 1);
            end
            dxi_dx = 1.0 / dx_dxi;

            for aa = 1 : n_en
                u_l    = u_l    + u_ele(aa) * PolyShape(pp, aa, xi(qua), 0);
                u_l_x  = u_l_x  + u_ele(aa) * PolyShape(pp, aa, xi(qua), 1) * dxi_dx;
            end
            
            y_l    = x_l^5;
            y_l_x  = 5*x_l^4;

            el2_n(ee) = el2_n(ee) + weight(qua) * ( u_l - y_l )^2 * dx_dxi;
            el2_d(ee) = el2_d(ee) + weight(qua) * ( y_l )^2 * dx_dxi;

            eh1_n(ee) = eh1_n(ee) + weight(qua) * ( u_l_x - y_l_x )^2 * dxi_dx;
            eh1_d(ee) = eh1_d(ee) + weight(qua) * ( y_l_x )^2 * dxi_dx;

        end

        % Assembly local integral
        El2_n = El2_n + el2_n(ee);
        El2_d = El2_d + el2_d(ee);
        Eh1_n = Eh1_n + eh1_n(ee);
        Eh1_d = Eh1_d + eh1_d(ee);

        for ll = 1 : n_sam_end
            x_l = 0.0;
            u_l = 0.0;
            u_l_x = 0.0;
            dx_dxi= 0.0;
            for aa = 1 : n_en
                dx_dxi = dx_dxi + x_ele(aa) * PolyShape(pp, aa, xi_sam(ll), 1);
            end
            dxi_dx = 1 / dx_dxi;
            for aa = 1 : n_en
                x_l   = x_l   + x_ele(aa) * PolyShape(pp, aa, xi_sam(ll), 0);
                u_l   = u_l   + u_ele(aa) * PolyShape(pp, aa, xi_sam(ll), 0);
                u_l_x = u_l_x + u_ele(aa) * PolyShape(pp, aa, xi_sam(ll), 1) * dxi_dx;
            end

            x_sam( (ee-1)*n_sam + ll ) = x_l;
            u_sam( (ee-1)*n_sam + ll ) = u_l;
            y_sam( (ee-1)*n_sam + ll ) = x_l^5;
            u_x_sam( (ee-1)*n_sam + ll ) = u_l_x;
            y_x_sam( (ee-1)*n_sam + ll ) = 5*x_l^4;
        end
    end

    El2(s) = sqrt(El2_n/El2_d);
    Eh1(s) = sqrt(Eh1_n/Eh1_d); %not true
end

log_h   = log(hh);
log_El2 = log(El2);
log_Eh1 = log(Eh1);

plot(log_h,log_El2,'b','LineWidth',2);
hold on
plot(log_h,log_Eh1,'r','LineWidth',2);

k_1=zeros(7,1); % slope of El2
k_2=zeros(7,1); % slope of Eh1

for s=1:7
    k_1(s) = (log_El2(s+1)-log_El2(s))/(log_h(s+1)-log_h(s));
    k_2(s) = (log_Eh1(s+1)-log_Eh1(s))/(log_h(s+1)-log_h(s));
end

K_1 = 0.0;
K_2 = 0.0;
for s=1:7
    K_1 = K_1 + k_1(s);
    K_2 = K_2 + k_2(s);
end
k_1
K_1=K_1/7
k_2
K_2=K_2/7

% plot(x_sam, y_sam, '-k','LineWidth',2);
% hold on;
% plot(x_sam, u_sam, '--r','LineWidth',2);
% 
% plot(x_sam, u_x_sam, '-b','LineWidth',2); % ？
% hold on;
% plot(x_sam, y_x_sam, '-k','LineWidth',2);

































% EOF