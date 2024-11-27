clear all; clc; clf;  % clean the memory, screen, and figure

% Problem definition
f = @(x) -20*x.^3; % f(x) is the source
g = 1.0;           % u    = g  at x = 1
h = 0.0;           % -u,x = h  at x = 0

% Setup the mesh
pp   = 1;              % polynomial degree
n_en = pp + 1;         % number of element or local nodes
n_el = 2;              % number of elements
n_np = n_el * pp + 1;  % number of nodal points
n_eq = n_np - 1;       % number of equations
n_int = 10;

hh = 1.0 / (n_np - 1); % space between two adjacent nodes
x_coor = 0 : hh : 1;   % nodal coordinates for equally spaced nodes

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

disp = [d_temp; g];

% Postprocessing: visualization

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

eh2_n = zeros(n_el,1); % numerator
eh2_d = el2_n;         % denominator
Eh2_n = 0;
Eh2_d = 0;

for ee = 1 : n_el
  x_ele = x_coor( IEN(ee, :) );
  u_ele = disp( IEN(ee, :) );

  if ee == n_el
    n_sam_end = n_sam+1;
  else
    n_sam_end = n_sam;
  end

  % quadrature loop
  for qua = 1 : n_int    
    dx_dxi = 0.0;
    x_l    = 0.0;
    u_l    = 0.0;
    u_l_x  = 0.0;
    for aa = 1 : n_en
      x_l    = x_l    + x_ele(aa) * PolyShape(pp, aa, xi(qua), 0); % x for quadrature rule
      dx_dxi = dx_dxi + x_ele(aa) * PolyShape(pp, aa, xi(qua), 1);
      u_l    = u_l    + u_ele(aa) * PolyShape(pp, aa, xi(qua), 0);
      u_l_x  = u_l_x  + u_ele(aa) * PolyShape(pp, aa, xi(qua), 1);
    end
    dxi_dx = 1.0 / dx_dxi;
    y_l    = x_l^5;
    y_l_x  = 5*x_l^4;

    el2_n(ee) = el2_n(ee) + weight(qua) * ( u_l - y_l )^2 * dx_dxi;
    el2_d(ee) = el2_d(ee) + weight(qua) * ( y_l )^2 * dx_dxi;

    eh2_n(ee) = eh2_n(ee) + weight(qua) * ( u_l_x - y_l_x )^2 * dxi_dx;
    eh2_d(ee) = eh2_d(ee) + weight(qua) * ( y_l_x )^2 * dxi_dx;

  end

  % Assembly local integral
  El2_n = El2_n + el2_n(ee); 
  El2_d = El2_d + el2_d(ee);
  Eh2_n = Eh2_n + eh2_n(ee);
  Eh2_d = Eh2_d + eh2_d(ee);

  for ll = 1 : n_sam_end
    x_l = 0.0;
    u_l = 0.0;
    u_l_x = 0.0;
    for aa = 1 : n_en
      x_l   = x_l   + x_ele(aa) * PolyShape(pp, aa, xi_sam(ll), 0);
      u_l   = u_l   + u_ele(aa) * PolyShape(pp, aa, xi_sam(ll), 0);
      u_l_x = u_l_x + u_ele(aa) * PolyShape(pp, aa, xi_sam(ll), 1);
    end

    x_sam( (ee-1)*n_sam + ll ) = x_l;
    u_sam( (ee-1)*n_sam + ll ) = u_l;
    y_sam( (ee-1)*n_sam + ll ) = x_l^5;
    u_x_sam( (ee-1)*n_sam + ll ) = u_l_x;
    y_x_sam( (ee-1)*n_sam + ll ) = 5*x_l^4;
  end
end

plot(x_sam, u_sam, '-r','LineWidth',2);
hold on;
plot(x_sam, y_sam, '-k','LineWidth',2);

plot(x_sam, u_x_sam, '-b','LineWidth',2); % ？
hold on;
plot(x_sam, y_x_sam, '-k','LineWidth',2);

El2 = sqrt(El2_n/El2_d);
Eh2 = sqrt(Eh2_n/Eh2_d);































% EOF