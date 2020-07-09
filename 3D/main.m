close all, clear all, clc
%{
This codes solves a given nonlinear 3D structure of a single material with
hyperelasticity using Newton-Raphson iteration. The avaliable elemnent is
the 20 node hexaedrum. The programm returns the final Cauchy tensor and
deformed shape.
                        19-----18-----17
                       /|            /|
                      / |           / |
                     /  12         /  11
                    20  |         16  |
                   /    |        /    |
                  /     7------6/-----5
                 /     /       /     /
8-node element: 13-----14-----15    /
                |    /        |    /
                |   8         |   4
                9  /          10 /
                | /           | /
                |/            |/
                1------2------3

MADE BY: Juan Sebastián Delgado Trujillo 
%}
X = 1; Y = 2; Z = 3;   % variables for the reading of the code
%% Problem parameters - INPUT:

E  = 210e6;    % Young modulus [kPa]
nu = 0.25;     % Poisson ratio [-]

% nodal coordinates ( x [m], y [m]):
initial = load('3D_coord.txt');

% nodal forces (magnitude [kN], node, axis):
%force = load('3D_force.txt');
%force = load('3D_force_tension.txt');
force = load('3D_force_torsion.txt');

n_steps = 50;    % number of force factors 

% supports (node, axis):
supports = load('3D_supports.txt');
        
% correspondence matrix (row = element, column = local nodal numeration):
LaG = load('3D_LAG.txt');
 
nod_ele = 20;       % nodes per element(quadratic = 20)

it_max  = 50;      % maximum number of iterations
tol     = 1e-5;    % tolerance
        
%% Problem's parameters - MAIN:

n_nod   = size(initial,  1);    % number of nodes
n_el    = size(LaG,      1);    % number of elements
n_dof   = 3*n_nod;              % number of degree of freedom (3 per node)
n_force = size(force,    1);    % number of nodal forces
n_restr = size(supports, 1);    % number of restricted nodes

% degree of fredomm matrix (row = node, column = DOF):
dof = [(1:3:n_dof)', (2:3:n_dof)', (3:3:n_dof)'];

% Lamme's constant:
lambda = (nu*E)/((1 + nu)*(1 - 2*nu));
miu    = E/(2*(1 + nu));

steps  = linspace(0, 1, n_steps);    % force factors

% restricted(cc) and (dd) DOFs are grouped in vectors:
cc = zeros(n_restr, 1);
for sp = 1:n_restr
    cc(sp) = dof(supports(sp, 1), supports(sp, 2));
end
dd = setdiff(1:n_dof, cc);

%% plot:
figure
hold on
fig_1 = [1:8, 1];
fig_2 = [13:20, 13];
fig_3 = [1, 2, 3, 10, 15, 14, 13,  9, 1];
fig_4 = [7, 6, 5, 11, 17, 18, 19, 12, 7];
for e = 1:n_el
    plot3(initial(LaG(e, fig_1), X), initial(LaG(e, fig_1), Y),...
          initial(LaG(e, fig_1), Z), 'b-o', 'LineWidth', 2)
    plot3(initial(LaG(e, fig_2), X), initial(LaG(e, fig_2), Y),...
          initial(LaG(e, fig_2), Z), 'b-o', 'LineWidth', 2)
    plot3(initial(LaG(e, fig_3), X), initial(LaG(e, fig_3), Y),...
          initial(LaG(e, fig_3), Z), 'b-o', 'LineWidth', 2)
    plot3(initial(LaG(e, fig_4), X), initial(LaG(e, fig_4), Y),...
          initial(LaG(e, fig_4), Z), 'b-o', 'LineWidth', 2)
     
     cg = mean(initial(LaG(e, :), :));
     text(cg(1), cg(2), cg(3), num2str(e))
end
grid on
axis equal
view([1, 1, 1])
title('Mesh')

%% FEM parameters:

% the interpolation functions, their derivatives w.r.t its arguments and
% the integration method (Gauss-Legendre quadratures) are defined:
[N, dN_ds, dN_dt, dN_dr, pg, w] = parameters(nod_ele);

n_gp = length(w);    % number of integration points

%% initial configuration:
% as the deformation matriz and the jacobian in the initial configuration 
% and at all Gauss points reamins constant through the iterations they are
% computed at once

B      = cell(n_el, n_gp);
J_0    = cell(n_el, n_gp);
stress = cell(n_el, n_gp);

% in each element:
for el = 1:n_el
    
    %in each Gauss point:
    for pt = 1:n_gp
        
        % 1) the Jacobian matrix is calculated
        j_aux = jacobian(initial(LaG(el, :), :), pg(pt, :));
        
        % 2) the derivatives of the interpolation functions w.r.t. the
        % initial coordinates are computed
        dNdx0 = (j_aux')\[dN_ds(pg(pt, 1), pg(pt, 2), pg(pt, 3))';
                          dN_dt(pg(pt, 1), pg(pt, 2), pg(pt, 3))';
                          dN_dr(pg(pt, 1), pg(pt, 2), pg(pt, 3))'];
                       
        % 3) the deformation matrix is made up
        B_aux = zeros(3*nod_ele, 9);
        for f = 1:nod_ele
            B_aux((3*f - 2):(3*f), :) =...
                [dNdx0(X, f), dNdx0(Y, f), dNdx0(Z, f),           0,           0,           0,           0,           0,           0;
                           0,           0,           0, dNdx0(X, f), dNdx0(Y, f), dNdx0(Z, f),           0,           0,           0;
                           0,           0,           0,           0,           0,           0, dNdx0(X, f), dNdx0(Y, f), dNdx0(Z, f)];
        end
                    
        % 4) the Jacobian matrix and the deformation matrix are stored
        J_0{el, pt} = j_aux;
        B{  el, pt} = B_aux;
    end
end

%% main:

% the current coordinates are initialized as the initial
current = initial;

for ld = 1:n_steps
    
    f_factor = steps(ld);    % the load factor for the current
    
    it   = 0;
    conv = 1;
    uu   = sparse(n_dof, 1);
    
    while (conv > tol) && (it < it_max)
        
        K    = sparse(n_dof, n_dof);
        Ri   = sparse(n_dof, 1);
        Re   = sparse(n_dof, 1);
        
        % the tangent stiffness matrix and the internal reaction vector are
        % assessed
            % in each element:
        for el = 1:n_el

            dof_el = reshape(dof(LaG(el, :), :)', 3*nod_ele, 1); % EF's DOF
            init_c = initial(LaG(el, :), :);     % EF's intial coordinates
            curr_c = current(LaG(el, :), :);     % EF's current coordinates
            
            k      = zeros(3*nod_ele);
            ri     = zeros(3*nod_ele, 1);

                % in each Gauss point:
            for pt = 1:n_gp

                % 1) the deformation gradient, deformation matrix,
                % constitutive matrix and the second PK tensor in the
                % current step are evaluated
                [F, F_b]    = def_grad(init_c, curr_c, pg(pt, :));
                B_pg        = B{el, pt};
                [C, S, S_b] = hyperelast(lambda, miu, F);

                % 2) the stiffness matrices (material and geometric) and
                % the internal reaction vector are calculated
                kc    =  B_pg*F_b'*C*F_b*B_pg'*det(J_0{el, pt})*w(pt);
                ks    =  B_pg*S_b*B_pg'       *det(J_0{el, pt})*w(pt);
                r_int = -B_pg*F_b'*S          *det(J_0{el, pt})*w(pt);

                % 3) the stiffness and the internal reactions are added to
                % the respective element array
                k  = k  + kc + ks;
                ri = ri + r_int;
            end

            % then, the stiffness matrix and the reaction vector of the
            % element are added in the global arrays
            K(dof_el, dof_el) = K(dof_el, dof_el) + k;
            Ri(dof_el)        = Ri(dof_el)        + ri;
        end

        % the external loads vector is ensambled
        for f = 1:n_force
            Re(dof(force(f, 2), force(f, 3))) = force(f, 1)*f_factor;
        end

        force_res = Re + Ri;    % force residual vector

        % the stiffness matrix and the force resisdual vector are
        % restricted by the supports
        K_dd  = K(dd, dd);
        fr_dd = force_res(dd);
        
        % the displacements are calculated by solving the algebraic system
        u = K_dd\fr_dd;
        uu(dd) = u;

        current = current + reshape(uu, 3, n_nod)';   % current coordinates

        % the internal reactions in the current configuration are assessed:
        Ri = sparse(n_dof, 1);
            % in each element:
        for el = 1:n_el

            dof_el = reshape(dof(LaG(el, :), :)', 3*nod_ele, 1); % EF's DOF
            init_c = initial(LaG(el, :), :);     % EF's intial coordinates
            curr_c = current(LaG(el, :), :);     % EF's current coordinates
            
            ri     = zeros(3*nod_ele, 1);

                % in each Gauss point:
            for pt = 1:n_gp

                % 1) the current Cauchy stresses is calculated
                [F, ~] = def_grad(init_c, curr_c, pg(pt, :));
                b      = F*F';
                Cauchy = log(det(F))*lambda/(det(F))*eye(3) + ...
                         miu/(det(F))*(b - eye(3));
                     
                sigma = [Cauchy(1, 1); Cauchy(2, 2); Cauchy(3, 3);
                         Cauchy(1, 2); Cauchy(2, 3); Cauchy(1, 3)];
                
                
                % 2) the deformation matrix in the current configuration is
                % computed
                j_aux = jacobian(current(LaG(el, :), :), pg(pt, :));
                dN_dx = (j_aux')\[dN_ds(pg(pt, 1), pg(pt, 2), pg(pt, 3))';
                                  dN_dt(pg(pt, 1), pg(pt, 2), pg(pt, 3))'; 
                                  dN_dr(pg(pt, 1), pg(pt, 2), pg(pt, 3))'];

                B_l = zeros(3*nod_ele, 6);
                for f = 1:nod_ele
                    B_l((3*f - 2):3*f, :) =...
                        [dN_dx(X, f),           0,           0, dN_dx(Y, f),           0, dN_dx(Z, f);
                                   0, dN_dx(Y, f),           0, dN_dx(X, f), dN_dx(Z, f),           0;
                                   0,           0, dN_dx(Z, f),           0, dN_dx(Y, f), dN_dx(X, f)];
                end

                % 3) the internal reaction vector are calculated and added
                % to the element vector
                ri = ri + B_l*sigma*det(j_aux)*w(pt);
                
                % 4) the Cauhcy stresses are saved if the current iteration
                % is the last
                if ld == n_steps
                    stress{el, pt} = sigma;
                end
            end

            % then, the reaction vector of the element are added to the
            % vector
            Ri(dof_el) = Ri(dof_el) + ri;
        end

        % finally, the convergence of the current solution is checked
        resi = Re - Ri;                              % residual
        conv = norm(resi(dd))/(1 + norm(Re(dd)));    % convergence
        
        it = it + 1;
        
    end
end

%% Report:

figure
hold on
fig_1 = [1:8, 1];
fig_2 = [13:20, 13];
fig_3 = [1, 2, 3, 10, 15, 14, 13,  9, 1];
fig_4 = [7, 6, 5, 11, 17, 18, 19, 12, 7];
for e = 1:n_el
    plot3(current(LaG(e, fig_1), X), current(LaG(e, fig_1), Y),...
          current(LaG(e, fig_1), Z), 'b-', 'LineWidth', 2)
    plot3(current(LaG(e, fig_2), X), current(LaG(e, fig_2), Y),...
          current(LaG(e, fig_2), Z), 'b-', 'LineWidth', 2)
    plot3(current(LaG(e, fig_3), X), current(LaG(e, fig_3), Y),...
          current(LaG(e, fig_3), Z), 'b-', 'LineWidth', 2)
    plot3(current(LaG(e, fig_4), X), current(LaG(e, fig_4), Y),...
          current(LaG(e, fig_4), Z), 'b-', 'LineWidth', 2)
      
      
    plot3(initial(LaG(e, fig_1), X), initial(LaG(e, fig_1), Y),...
          initial(LaG(e, fig_1), Z), 'k--', 'LineWidth', 1)
    plot3(initial(LaG(e, fig_2), X), initial(LaG(e, fig_2), Y),...
          initial(LaG(e, fig_2), Z), 'k--', 'LineWidth', 1)
    plot3(initial(LaG(e, fig_3), X), initial(LaG(e, fig_3), Y),...
          initial(LaG(e, fig_3), Z), 'k--', 'LineWidth', 1)
    plot3(initial(LaG(e, fig_4), X), initial(LaG(e, fig_4), Y),...
          initial(LaG(e, fig_4), Z), 'k--', 'LineWidth', 1)
end
grid on
axis equal
title('Deformed mesh')
view([1, 1, 1])

to_gid(initial, current, LaG, stress, nod_ele)    % resultados en GID