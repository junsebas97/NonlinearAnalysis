function [F, F_bar] = def_grad(init_coord, curr_coord, GP)
%{
This function calcules the deformation gradient between the initial and the
current configartion evaluated in a Gauss point of the element

Args:
init_coord: initial coordinates of the FE [x0, y0]
curr_coord: coordinates of the current configuration of the FE [xc, yc]
pg: Gauss point

Outputs: the deformation gradient matrix
%}

s = GP(1);        t = GP(2);

% derivatives of the interpolation functions:
n_nod = size(init_coord, 1);
if n_nod == 4
    dN_ds = [(t - 1)/4;  (1 - t)/4; (t + 1)/4; -(t + 1)/4];
    dN_dt = [(s - 1)/4; -(s + 1)/4; (s + 1)/4;  (1 - s)/4];
elseif n_nod == 8
    dN_ds = [-((2*s + t)*(t - 1))/4;           s*(t - 1);
             -((2*s - t)*(t - 1))/4;         1/2 - t^2/2;
              ((2*s + t)*(t + 1))/4;          -s*(t + 1);
              ((2*s - t)*(t + 1))/4;         t^2/2 - 1/2];

    dN_dt = [-((s + 2*t)*(s - 1))/4;         s^2/2 - 1/2;
             -((s - 2*t)*(s + 1))/4;          -t*(s + 1);
              ((s + 2*t)*(s + 1))/4;         1/2 - s^2/2;
              ((s - 2*t)*(s - 1))/4;           t*(s - 1)];
end

% the Jacobian matrix in the inital and current configuration are
% calculated at the given point
dx0_ds = init_coord(:, 1)'*dN_ds;        dx0_dt = init_coord(:, 1)'*dN_dt;
dy0_ds = init_coord(:, 2)'*dN_ds;        dy0_dt = init_coord(:, 2)'*dN_dt;

dx_ds = curr_coord(:, 1)'*dN_ds;         dx_dt  = curr_coord(:, 1)'*dN_dt;
dy_ds = curr_coord(:, 2)'*dN_ds;         dy_dt  = curr_coord(:, 2)'*dN_dt;

J0 = [dx0_ds, dx0_dt; dy0_ds, dy0_dt];
J  = [ dx_ds,  dx_dt;  dy_ds,  dy_dt];

% then, the deformation gradient tensor is assessed
F = (J0'\J')';

% Finally, it is writen in matrix notation
F_bar = [F(1, 1),       0, F(2, 1),       0;
               0, F(1, 2),       0, F(2, 2);
         F(1, 2), F(1, 1), F(2, 2), F(2, 1)];
end