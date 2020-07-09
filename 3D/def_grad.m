function [F, F_bar] = def_grad(i_coor, c_coor, GP)
%{
This function calcules the deformation gradient between the initial and the
current configartion evaluated in a Gauss point of the element

Args:
i_coor: initial coordinates of the FE [x0, y0]
c_coor: coordinates of the current configuration of the FE [xc, yc]
pg:     Gauss point

Returns: the deformation gradient tensor and matrix
%}
s = GP(1);        t = GP(2);        r = GP(3);

% derivatives of the interpolation functions:
n_nod = size(i_coor, 1);
[~, dNds, dNdt, dNdr, ~, ~] = parameters(n_nod);

% the Jacobian matrix in the inital and current configuration are
% calculated at the given point
dx0ds = i_coor(:, 1)'*dNds(s, t, r);    dxds = c_coor(:, 1)'*dNds(s, t, r);
dx0dt = i_coor(:, 1)'*dNdt(s, t, r);    dxdt = c_coor(:, 1)'*dNdt(s, t, r);
dx0dr = i_coor(:, 1)'*dNdr(s, t, r);    dxdr = c_coor(:, 1)'*dNdr(s, t, r);

dy0ds = i_coor(:, 2)'*dNds(s, t, r);    dyds = c_coor(:, 2)'*dNds(s, t, r);
dy0dt = i_coor(:, 2)'*dNdt(s, t, r);    dydt = c_coor(:, 2)'*dNdt(s, t, r);
dy0dr = i_coor(:, 2)'*dNdr(s, t, r);    dydr = c_coor(:, 2)'*dNdr(s, t, r);

dz0ds = i_coor(:, 3)'*dNds(s, t, r);    dzds = c_coor(:, 3)'*dNds(s, t, r);
dz0dt = i_coor(:, 3)'*dNdt(s, t, r);    dzdt = c_coor(:, 3)'*dNdt(s, t, r);
dz0dr = i_coor(:, 3)'*dNdr(s, t, r);    dzdr = c_coor(:, 3)'*dNdr(s, t, r);


J0 = [dx0ds, dx0dt, dx0dr;
      dy0ds, dy0dt, dy0dr;
      dz0ds, dz0dt, dz0dr];
  
J  = [ dxds,  dxdt,  dxdr;
       dyds,  dydt,  dydr;
       dzds,  dzdt,  dzdr];

% then, the deformation gradient tensor is assessed
F = (J0'\J')';

% Finally, it is writen in matrix notation
F_bar = [F(1, 1),       0,       0, F(2, 1),       0,       0, F(3, 1),       0,       0;
               0, F(1, 2),       0,       0, F(2, 2),       0,       0, F(3, 2),       0;
               0,       0, F(1, 3),       0,       0, F(2, 3),       0,       0, F(3, 3);
         F(1, 2), F(1, 1),       0, F(2, 2), F(2, 1),       0, F(3, 2), F(3, 1),       0;
               0, F(1, 3), F(1, 2),       0, F(2, 3), F(2, 2),       0, F(3, 3), F(3, 2);
         F(1, 3),       0, F(1, 1), F(2, 3),       0, F(2, 1), F(3, 3),       0, F(3, 1)];
end