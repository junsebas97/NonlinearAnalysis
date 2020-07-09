close all, clear all, clc
%{
This codes solves a given nonlinear structure of a single material with
plane stress analysis and hyperelasticity using Newton-Raphson iteration.
2 elemnents are avalaible 4-node bilinear or 8-node quadratic
The programm returns the final Cauchy tensor and deformed shape.

4-node element:     4--------3            8-node element: 7---6---5
                    |        |                            |       |
                    |        |                            8       4
                    |        |                            |       |
                    1--------2                            1---2---3

MADE BY: Juan Sebastián Delgado Trujillo 
%}
X = 1; Y = 2;    % variables for the reading of the code
%% Problem parameters - INPUT:

E  = 210e6;    % Young modulus [kPa]
nu = 0.25;     % Poisson ratio [-]
t  = 0.1;      % thickness     [m]

% nodal coordinates ( x [m], y [m]):
%initial = load('2D4N_coord.txt');
initial = load('2D8N_coord.txt');

% nodal forces (magnitude [kN], node, axis):
%force = load('2D4N_force.txt');
force = load('2D8N_force.txt');

n_steps = 20;    % number of force factors 

% supports (node, axis):
%supports = load('2D4N_supports.txt');
supports = load('2D8N_supports.txt');
        
% correspondence matrix (row = element, column = local nodal numeration):
%LaG = load('2D4N_LAG.txt');
LaG = load('2D8N_LAG.txt');
 
%nod_ele = 4;       % nodes per element(bi-linear = 4; quadratic = 8)
nod_ele = 8;

it_max  = 50;      % maximum number of iterations
tol     = 1e-5;    % tolerance
        
%% Problem's parameters - MAIN:

n_nod   = size(initial,  1);    % number of nodes
n_el    = size(LaG,      1);    % number of elements
n_dof   = 2*n_nod;              % number of degree of freedom (2 per node)
n_force = size(force,    1);    % number of nodal forces
n_restr = size(supports, 1);    % number of restricted nodes

% degree of fredomm matrix (row = node, column = DOF):
dof = [(1:2:n_dof)', (2:2:n_dof)'];

% Lamme's constant:
lambda = (nu*E)/((1 + nu)*(1 - 2*nu));
miu    = E/(2*(1 + nu));
gamma  = (2*miu)/(lambda + 2*miu);

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
for e = 1:n_el
    plot([initial(LaG(e, :), X); initial(LaG(e, 1), X)], ...
         [initial(LaG(e, :), Y); initial(LaG(e, 1), Y)], ...
         'b-o', 'LineWidth', 2)
     
     cg = mean(initial(LaG(e, :), :));
     text(cg(1), cg(2), num2str(e), 'FontSize', 12)
end
grid on
axis equal
title('Mesh')

%% FEM parameters:

% the interpolation functions, their derivatives w.r.t its arguments and
% the integration method (Gauss-Legendre quadratures) are defined:
[N, dN_ds, dN_dt, pg, w] = parameters(nod_ele);

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
        dNdx0 = (j_aux')\[dN_ds(pg(pt, 1), pg(pt, 2))';
                           dN_dt(pg(pt, 1), pg(pt, 2))'];
                       
        % 3) the deformation matrix is made up
        B_aux = zeros(2*nod_ele, 4);
        for f = 1:nod_ele
            B_aux((2*f - 1):(2*f), :) = [dNdx0(X, f), dNdx0(Y, f),            0,            0;
                                                   0,            0, dNdx0(X, f), dNdx0(Y, f)];
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

            dof_el = reshape(dof(LaG(el, :), :)', 2*nod_ele, 1); % EF's DOF
            init_c = initial(LaG(el, :), :);     % EF's intial coordinates
            curr_c = current(LaG(el, :), :);     % EF's current coordinates
            
            k      = zeros(2*nod_ele);
            ri     = zeros(2*nod_ele, 1);

                % in each Gauss point:
            for pt = 1:n_gp

                % 1) the deformation gradient, deformation matrix,
                % constitutive matrix and the second PK tensor in the
                % current step are evaluated
                [F, F_b]    = def_grad(init_c, curr_c, pg(pt, :));
                B_pg        = B{el, pt};
                [C, S, S_b] = hyperelast(lambda, miu, gamma, F);

                % 2) the stiffness matrices (material and geometric) and
                % the internal reaction vector are calculated
                kc    =  B_pg*F_b'*C*(B_pg*F_b')'*det(J_0{el, pt})*t*w(pt);
                ks    =  B_pg*S_b*B_pg'          *det(J_0{el, pt})*t*w(pt);
                r_int = -B_pg*F_b'*S             *det(J_0{el, pt})*t*w(pt);

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

        current = current + reshape(uu, 2, n_nod)';   % current coordinates

        % the internal reactions in the current configuration are assessed:
        Ri = sparse(n_dof, 1);
            % in each element:
        for el = 1:n_el

            dof_el = reshape(dof(LaG(el, :), :)', 2*nod_ele, 1); % EF's DOF
            init_c = initial(LaG(el, :), :);     % EF's intial coordinates
            curr_c = current(LaG(el, :), :);     % EF's current coordinates
            
            ri     = zeros(2*nod_ele, 1);

                % in each Gauss point:
            for pt = 1:n_gp

                % 1) the current Cauchy stresses and the current thickness
                % are calculated
                [F, ~] = def_grad(init_c, curr_c, pg(pt, :));
                b      = F*F';
                Cauchy = (log(det(F)^gamma))*lambda/(det(F)^gamma)*eye(2) + ...
                          miu/(det(F)^gamma)*(b - eye(2));

                t_c   = t*det(F)^(gamma - 1);
                sigma = [Cauchy(1, 1); Cauchy(2, 2); Cauchy(1, 2)];
                
                
                % 2) the deformation matrix in the current configuration is
                % computed
                j_aux = jacobian(current(LaG(el, :), :), pg(pt, :));
                dN_dx = (j_aux')\[dN_ds(pg(pt, 1), pg(pt, 2))';
                                  dN_dt(pg(pt, 1), pg(pt, 2))'];
                
                B_l = zeros(2*nod_ele, 3);
                for f = 1:nod_ele
                    B_l((2*f - 1):2*f, :) = [dN_dx(X, f),           0,  dN_dx(Y, f);
                                                       0, dN_dx(Y, f),  dN_dx(X, f)];
                end

                % 3) the internal reaction vector are calculated and added
                % to the element vector
                ri = ri + B_l*sigma*det(j_aux)*t_c*w(pt);
                
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
for e = 1:n_el
    plot([current(LaG(e, :), X); current(LaG(e, 1), X)], ...
         [current(LaG(e, :), Y); current(LaG(e, 1), Y)], ...
         'b-', 'LineWidth', 2)
     
    plot([initial(LaG(e, :), X); initial(LaG(e, 1), X)], ...
         [initial(LaG(e, :), Y); initial(LaG(e, 1), Y)], ...
         'k--', 'LineWidth', 1)
end
grid on
axis equal
title('Deformed mesh')

to_gid(initial, current, LaG, stress, nod_ele) % resultados en GID